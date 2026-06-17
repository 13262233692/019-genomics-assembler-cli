package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/debug"
	"sort"
	"sync"
	"sync/atomic"
	"time"

	"github.com/genomics-assembler/cli/assembler"
	"github.com/genomics-assembler/cli/debruijn"
	"github.com/genomics-assembler/cli/fasta"
	"github.com/genomics-assembler/cli/fastq"
)

var (
	inputFile  string
	outputFile string
	kmerSize   int
	minContig  int
	numWorkers int
	verbose    bool
	maxContig  int
	maxNodes   uint64
	maxEdges   uint64
)

func init() {
	flag.StringVar(&inputFile, "input", "-", "Input FASTQ file path (use - for stdin). Supports .gz files.")
	flag.StringVar(&inputFile, "i", "-", "Short for -input")
	flag.StringVar(&outputFile, "output", "contigs.fasta", "Output FASTA file path")
	flag.StringVar(&outputFile, "o", "contigs.fasta", "Short for -output")
	flag.IntVar(&kmerSize, "kmer", 31, "K-mer size for De Bruijn graph")
	flag.IntVar(&kmerSize, "k", 31, "Short for -kmer")
	flag.IntVar(&minContig, "min-length", 100, "Minimum contig length to output")
	flag.IntVar(&minContig, "m", 100, "Short for -min-length")
	flag.IntVar(&numWorkers, "workers", runtime.NumCPU(), "Number of worker goroutines")
	flag.IntVar(&numWorkers, "w", runtime.NumCPU(), "Short for -workers")
	flag.BoolVar(&verbose, "verbose", false, "Verbose output")
	flag.BoolVar(&verbose, "v", false, "Short for -verbose")
	flag.IntVar(&maxContig, "max-contig", 10_000_000, "Maximum contig length (deadlock protection)")
	flag.Uint64Var(&maxNodes, "max-nodes", 500_000_000, "Maximum graph node count (OOM protection)")
	flag.Uint64Var(&maxEdges, "max-edges", 2_000_000_000, "Maximum graph edge count (OOM protection)")
}

func main() {
	flag.Parse()

	if kmerSize < 3 {
		log.Fatal("k-mer size must be at least 3")
	}

	log.Printf("Genome Assembler CLI starting...")
	log.Printf("K-mer size: %d, Min contig: %d, Max contig: %d, Workers: %d",
		kmerSize, minContig, maxContig, numWorkers)
	log.Printf("Max nodes: %d, Max edges: %d", maxNodes, maxEdges)
	log.Printf("Input: %s, Output: %s", inputFile, outputFile)

	start := time.Now()
	var memStats runtime.MemStats
	readMemStats(&memStats)
	log.Printf("Initial memory: %.2f MB", float64(memStats.Alloc)/1024/1024)

	reader, err := fastq.OpenReader(inputFile)
	if err != nil {
		log.Fatalf("Failed to open input: %v", err)
	}
	if reader != os.Stdin {
		defer reader.Close()
	}

	graph := debruijn.NewGraphWithLimits(kmerSize, maxNodes, maxEdges)

	readChan := make(chan *fastq.Read, numWorkers*10)
	var wg sync.WaitGroup
	var fatalErr atomic.Value
	var addErrCount uint64

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for read := range readChan {
				if err := graph.AddSequence(read.Sequence); err != nil {
					atomic.AddUint64(&addErrCount, 1)
					fatalErr.Store(err)
					if verbose {
						log.Printf("AddSequence error: %v", err)
					}
				}
			}
		}()
	}

	parser := fastq.NewParser(reader)
	readCount := 0
	reportEvery := uint64(1_000_000)

	for {
		read, err := parser.Next()
		if err == io.EOF || read == nil {
			break
		}
		if err != nil {
			log.Fatalf("Error reading FASTQ: %v", err)
		}
		readChan <- read
		readCount++
		if verbose && uint64(readCount)%reportEvery == 0 {
			readMemStats(&memStats)
			log.Printf("Processed %d reads | nodes=%d edges=%d | mem=%.2f MB",
				readCount, graph.NodeCount(), graph.EdgeCount(),
				float64(memStats.Alloc)/1024/1024)
		}
		if fatalErr.Load() != nil {
			break
		}
	}
	close(readChan)
	wg.Wait()

	if errVal := fatalErr.Load(); errVal != nil {
		log.Printf("WARNING: graph construction encountered error: %v (skipped %d sequences)",
			errVal, atomic.LoadUint64(&addErrCount))
	}

	graphBuildTime := time.Since(start)
	log.Printf("Reads processed: %d", readCount)
	log.Printf("Graph: nodes=%d edges=%d", graph.NodeCount(), graph.EdgeCount())
	log.Printf("Graph built in %s", graphBuildTime)

	if verbose {
		log.Println("Verifying graph integrity...")
		collisions := graph.VerifyIntegrity()
		if collisions > 0 {
			log.Printf("INTEGRITY: %d collision(s) detected in double-hash signatures", collisions)
		} else {
			log.Println("INTEGRITY: double-hash signatures OK, zero collisions")
		}
	}

	readMemStats(&memStats)
	log.Printf("Pre-assembly memory: %.2f MB", float64(memStats.Alloc)/1024/1024)

	log.Println("Assembling contigs...")
	asm := assembler.New()
	asm.MaxContigLength = maxContig
	asm.Verbose = verbose
	contigs := asm.Assemble(graph, minContig)
	assembleTime := time.Since(start) - graphBuildTime
	log.Printf("Contigs assembled: %d (in %s)", len(contigs), assembleTime)

	totalBases := 0
	maxLen := 0
	n50 := 0
	lengths := make([]int, 0, len(contigs))
	for _, c := range contigs {
		totalBases += c.Length
		lengths = append(lengths, c.Length)
		if c.Length > maxLen {
			maxLen = c.Length
		}
	}
	n50 = computeN50(lengths, totalBases)
	log.Printf("Assembly stats: total_bases=%d max_len=%d N50=%d", totalBases, maxLen, n50)

	log.Printf("Writing contigs to %s...", outputFile)
	writer, err := fasta.NewFileWriter(outputFile)
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	for _, c := range contigs {
		seq := &fasta.Sequence{
			Header:   fmt.Sprintf("contig_%d len=%d", c.ID, c.Length),
			Sequence: c.Sequence,
		}
		if err := writer.Write(seq); err != nil {
			log.Fatalf("Failed to write contig: %v", err)
		}
	}
	if err := writer.Flush(); err != nil {
		log.Fatalf("Failed to flush output: %v", err)
	}

	readMemStats(&memStats)
	log.Printf("Final memory: %.2f MB", float64(memStats.Alloc)/1024/1024)

	debug.FreeOSMemory()
	totalTime := time.Since(start)
	log.Printf("Done! Total time: %s", totalTime)
}

func readMemStats(m *runtime.MemStats) {
	runtime.ReadMemStats(m)
}

func computeN50(lengths []int, totalBases int) int {
	if len(lengths) == 0 || totalBases == 0 {
		return 0
	}
	sort.Ints(lengths)
	target := totalBases / 2
	sum := 0
	for i := len(lengths) - 1; i >= 0; i-- {
		sum += lengths[i]
		if sum >= target {
			return lengths[i]
		}
	}
	return lengths[len(lengths)-1]
}
