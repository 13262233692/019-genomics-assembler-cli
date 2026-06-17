package main

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"sync"
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
}

func main() {
	flag.Parse()

	if kmerSize < 3 {
		log.Fatal("k-mer size must be at least 3")
	}

	log.Printf("Genome Assembler CLI starting...")
	log.Printf("K-mer size: %d, Min contig length: %d, Workers: %d", kmerSize, minContig, numWorkers)
	log.Printf("Input: %s, Output: %s", inputFile, outputFile)

	start := time.Now()

	reader, err := fastq.OpenReader(inputFile)
	if err != nil {
		log.Fatalf("Failed to open input: %v", err)
	}
	if reader != os.Stdin {
		defer reader.Close()
	}

	graph := debruijn.NewGraph(kmerSize)

	readChan := make(chan *fastq.Read, numWorkers*10)
	var wg sync.WaitGroup

	for i := 0; i < numWorkers; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			for read := range readChan {
				graph.AddSequence(read.Sequence)
			}
		}()
	}

	parser := fastq.NewParser(reader)
	readCount := 0
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
		if verbose && readCount%1000000 == 0 {
			log.Printf("Processed %d reads, %d nodes so far...", readCount, graph.NodeCount())
		}
	}
	close(readChan)
	wg.Wait()

	graphBuildTime := time.Since(start)
	log.Printf("Reads processed: %d", readCount)
	log.Printf("Graph nodes: %d", graph.NodeCount())
	log.Printf("Graph built in %s", graphBuildTime)

	log.Println("Assembling contigs...")
	contigs := assembler.Assemble(graph, minContig)
	assembleTime := time.Since(start) - graphBuildTime
	log.Printf("Contigs assembled: %d (in %s)", len(contigs), assembleTime)

	totalBases := 0
	maxLen := 0
	for _, c := range contigs {
		totalBases += c.Length
		if c.Length > maxLen {
			maxLen = c.Length
		}
	}
	log.Printf("Total bases in contigs: %d", totalBases)
	log.Printf("Max contig length: %d", maxLen)

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

	totalTime := time.Since(start)
	log.Printf("Done! Total time: %s", totalTime)
}
