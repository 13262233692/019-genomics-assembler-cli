package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime"
	"runtime/debug"
	"sort"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"time"

	"github.com/genomics-assembler/cli/assembler"
	"github.com/genomics-assembler/cli/cmap"
	"github.com/genomics-assembler/cli/debruijn"
	"github.com/genomics-assembler/cli/fasta"
	"github.com/genomics-assembler/cli/fastq"
	"github.com/genomics-assembler/cli/restriction"
	"github.com/genomics-assembler/cli/scaffold"
)

const usageText = `Genome Assembler CLI - High-performance NGS genome assembly tool

Usage:
  genomics-assembler <command> [options]

Commands:
  assemble    Assemble reads into contigs (De Bruijn graph)
  scaffold    Order contigs into scaffolds using optical mapping (CMAP)
  help        Show this help message

Use "genomics-assembler <command> -h" for more information about a command.
`

func main() {
	if len(os.Args) < 2 {
		fmt.Print(usageText)
		os.Exit(1)
	}

	cmd := os.Args[1]
	os.Args = append([]string{os.Args[0]}, os.Args[2:]...)

	switch cmd {
	case "assemble":
		runAssemble()
	case "scaffold":
		runScaffold()
	case "help", "-h", "--help":
		fmt.Print(usageText)
	default:
		fmt.Printf("Unknown command: %s\n\n", cmd)
		fmt.Print(usageText)
		os.Exit(1)
	}
}

func runAssemble() {
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

	fs := flag.NewFlagSet("assemble", flag.ExitOnError)
	fs.StringVar(&inputFile, "input", "-", "Input FASTQ file path (use - for stdin). Supports .gz files.")
	fs.StringVar(&inputFile, "i", "-", "Short for -input")
	fs.StringVar(&outputFile, "output", "contigs.fasta", "Output FASTA file path")
	fs.StringVar(&outputFile, "o", "contigs.fasta", "Short for -output")
	fs.IntVar(&kmerSize, "kmer", 31, "K-mer size for De Bruijn graph")
	fs.IntVar(&kmerSize, "k", 31, "Short for -kmer")
	fs.IntVar(&minContig, "min-length", 100, "Minimum contig length to output")
	fs.IntVar(&minContig, "m", 100, "Short for -min-length")
	fs.IntVar(&numWorkers, "workers", runtime.NumCPU(), "Number of worker goroutines")
	fs.IntVar(&numWorkers, "w", runtime.NumCPU(), "Short for -workers")
	fs.BoolVar(&verbose, "verbose", false, "Verbose output")
	fs.BoolVar(&verbose, "v", false, "Short for -verbose")
	fs.IntVar(&maxContig, "max-contig", 10_000_000, "Maximum contig length (deadlock protection)")
	fs.Uint64Var(&maxNodes, "max-nodes", 500_000_000, "Maximum graph node count (OOM protection)")
	fs.Uint64Var(&maxEdges, "max-edges", 2_000_000_000, "Maximum graph edge count (OOM protection)")
	fs.Parse(os.Args[1:])

	if kmerSize < 3 {
		log.Fatal("k-mer size must be at least 3")
	}

	log.Printf("Genome Assembler CLI - assemble command starting...")
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

func runScaffold() {
	var (
		contigFile   string
		cmapFile     string
		outputFile   string
		enzymeName   string
		tolerance    float64
		minMatches   int
		minScore     float64
		gapSize      int
		verbose      bool
	)

	fs := flag.NewFlagSet("scaffold", flag.ExitOnError)
	fs.StringVar(&contigFile, "contigs", "", "Input contigs FASTA file (required)")
	fs.StringVar(&contigFile, "c", "", "Short for -contigs")
	fs.StringVar(&cmapFile, "cmap", "", "Input optical map CMAP file (required)")
	fs.StringVar(&cmapFile, "m", "", "Short for -cmap")
	fs.StringVar(&outputFile, "output", "scaffolds.fasta", "Output scaffolds FASTA file")
	fs.StringVar(&outputFile, "o", "scaffolds.fasta", "Short for -output")
	fs.StringVar(&enzymeName, "enzyme", "DLE-1", "Restriction enzyme name (DLE-1, BspQI, BbvCI, BssSI, SwaI, PacI)")
	fs.StringVar(&enzymeName, "e", "DLE-1", "Short for -enzyme")
	fs.Float64Var(&tolerance, "tolerance", scaffold.DefaultTolerance, "Site position alignment tolerance (bp)")
	fs.IntVar(&minMatches, "min-matches", scaffold.DefaultMinMatches, "Minimum matched sites for alignment")
	fs.Float64Var(&minScore, "min-score", scaffold.DefaultMinScore, "Minimum normalized alignment score")
	fs.IntVar(&gapSize, "gap-size", scaffold.DefaultGapSize, "Default gap size (N's) between contigs")
	fs.BoolVar(&verbose, "verbose", false, "Verbose output")
	fs.BoolVar(&verbose, "v", false, "Short for -verbose")
	fs.Parse(os.Args[1:])

	if contigFile == "" || cmapFile == "" {
		fmt.Println("Error: -contigs and -cmap are required")
		fs.PrintDefaults()
		os.Exit(1)
	}

	enzyme, ok := restriction.FindEnzyme(enzymeName)
	if !ok {
		log.Fatalf("Unknown enzyme: %s. Supported: DLE-1, BspQI, BbvCI, BssSI, SwaI, PacI", enzymeName)
	}

	log.Printf("Genome Assembler CLI - scaffold command starting...")
	log.Printf("Contigs: %s, CMAP: %s, Enzyme: %s", contigFile, cmapFile, enzyme.Name)
	log.Printf("Tolerance: %.0f bp, Min matches: %d, Min score: %.2f", tolerance, minMatches, minScore)
	log.Printf("Output: %s", outputFile)

	start := time.Now()
	var memStats runtime.MemStats
	readMemStats(&memStats)
	log.Printf("Initial memory: %.2f MB", float64(memStats.Alloc)/1024/1024)

	log.Println("Loading CMAP...")
	cmapData, err := cmap.OpenCMAP(cmapFile)
	if err != nil {
		log.Fatalf("Failed to open CMAP: %v", err)
	}
	log.Printf("CMAP loaded: %d contigs, %d total sites",
		cmapData.ContigCount(), cmapData.TotalSites())

	log.Println("Digesting contigs...")
	contigSeqs, err := loadFastaSequences(contigFile)
	if err != nil {
		log.Fatalf("Failed to load contigs: %v", err)
	}
	log.Printf("Loaded %d contigs from FASTA", len(contigSeqs))

	digests := make([]*restriction.ContigDigest, 0, len(contigSeqs))
	totalSites := 0
	for i, seq := range contigSeqs {
		d := restriction.DigestSequence(i, seq, enzyme)
		digests = append(digests, d)
		totalSites += d.SiteCount()
	}
	log.Printf("Digested %d contigs, %d total restriction sites", len(digests), totalSites)

	log.Println("Aligning contigs to optical map...")
	scaffolder := scaffold.NewScaffolder()
	scaffolder.Tolerance = tolerance
	scaffolder.MinMatches = minMatches
	scaffolder.MinScore = minScore
	scaffolder.GapSize = gapSize
	scaffolder.Verbose = verbose

	readMemStats(&memStats)
	log.Printf("Pre-alignment memory: %.2f MB", float64(memStats.Alloc)/1024/1024)

	scaffolds := scaffolder.Run(digests, cmapData)
	scaffoldTime := time.Since(start)
	log.Printf("Scaffolding complete in %s: %d scaffolds built", scaffoldTime, len(scaffolds))

	totalBases := 0
	maxLen := 0
	n50 := 0
	totalGaps := 0
	lengths := make([]int, 0, len(scaffolds))
	for _, sc := range scaffolds {
		totalBases += sc.Length
		totalGaps += sc.TotalGapBases
		lengths = append(lengths, sc.Length)
		if sc.Length > maxLen {
			maxLen = sc.Length
		}
	}
	n50 = computeN50(lengths, totalBases)
	log.Printf("Scaffold stats: total_bases=%d max_len=%d N50=%d gap_bases=%d (%.1f%%)",
		totalBases, maxLen, n50, totalGaps, 100.0*float64(totalGaps)/float64(totalBases))

	log.Printf("Writing scaffolds to %s...", outputFile)
	writer, err := fasta.NewFileWriter(outputFile)
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	for _, sc := range scaffolds {
		header := fmt.Sprintf("scaffold_%d len=%d contigs=%d map=%s",
			sc.ID, sc.Length, sc.ContigCount, sc.MapID)
		seq := &fasta.Sequence{
			Header:   header,
			Sequence: sc.Sequence,
		}
		if err := writer.Write(seq); err != nil {
			log.Fatalf("Failed to write scaffold: %v", err)
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

func loadFastaSequences(path string) ([]string, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	seqs := make([]string, 0)
	var currentSeq strings.Builder
	inSeq := false

	scanner := bufio.NewScanner(f)
	scanner.Buffer(make([]byte, 1024*1024), 100*1024*1024)

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		if strings.HasPrefix(line, ">") {
			if inSeq {
				seqs = append(seqs, currentSeq.String())
				currentSeq.Reset()
			}
			inSeq = true
		} else if inSeq {
			currentSeq.WriteString(line)
		}
	}
	if inSeq && currentSeq.Len() > 0 {
		seqs = append(seqs, currentSeq.String())
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return seqs, nil
}

func computeScaffoldN50(lengths []int, totalBases int) int {
	return computeN50(lengths, totalBases)
}

func mustParseInt(s string) int {
	v, _ := strconv.Atoi(s)
	return v
}
