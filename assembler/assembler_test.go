package assembler

import (
	"testing"
	"time"

	"github.com/genomics-assembler/cli/debruijn"
)

func TestAssembleSimple(t *testing.T) {
	g := debruijn.NewGraph(4)
	_ = g.AddSequence("ACGTGACTAGCT")

	contigs := Assemble(g, 5)
	if len(contigs) == 0 {
		t.Fatal("Expected at least one contig")
	}

	totalBases := 0
	for _, c := range contigs {
		totalBases += c.Length
		t.Logf("Contig %d: len=%d, seq=%s", c.ID, c.Length, c.Sequence)
	}
	if totalBases < 10 {
		t.Errorf("Expected at least 10 total bases, got %d", totalBases)
	}
}

func TestAssembleMinLength(t *testing.T) {
	g := debruijn.NewGraph(4)
	_ = g.AddSequence("ACGT")

	contigs := Assemble(g, 100)
	if len(contigs) != 0 {
		t.Errorf("Expected 0 contigs with min length 100, got %d", len(contigs))
	}

	contigs = Assemble(g, 2)
	if len(contigs) == 0 {
		t.Error("Expected at least 1 contig with min length 2")
	}
}

func TestAssembleMultipleContigs(t *testing.T) {
	g := debruijn.NewGraph(4)
	_ = g.AddSequence("AAGCTTGC")
	_ = g.AddSequence("CCGGAATT")

	contigs := Assemble(g, 5)
	if len(contigs) < 2 {
		for i, c := range contigs {
			t.Logf("Contig %d: len=%d, seq=%s", i, c.Length, c.Sequence)
		}
		t.Errorf("Expected at least 2 contigs, got %d", len(contigs))
	}
}

func TestAssembleRepetitiveSequenceDoesNotHang(t *testing.T) {
	done := make(chan bool, 1)
	go func() {
		g := debruijn.NewGraph(5)
		rep := "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
		_ = g.AddSequence(rep)
		asm := New()
		asm.MaxContigLength = 1000
		asm.MaxStepsPerWalk = 10000
		contigs := asm.Assemble(g, 5)
		t.Logf("Repetitive sequence produced %d contigs", len(contigs))
		for _, c := range contigs {
			t.Logf("  Contig %d: len=%d", c.ID, c.Length)
		}
		done <- true
	}()

	select {
	case <-done:
	case <-time.After(5 * time.Second):
		t.Fatal("Assembler hung on repetitive sequence (deadlock!)")
	}
}

func TestAssemblerDeadlockProtection(t *testing.T) {
	done := make(chan bool, 1)
	go func() {
		g := debruijn.NewGraph(5)
		_ = g.AddSequence("AAAAAAAAAA")
		_ = g.AddSequence("TTTTTTTTTT")
		_ = g.AddSequence("CCCCCCCCCC")
		_ = g.AddSequence("GGGGGGGGGG")

		asm := New()
		asm.MaxContigLength = 500
		asm.MaxStepsPerWalk = 1000
		contigs := asm.Assemble(g, 2)
		t.Logf("Max-stress test produced %d contigs", len(contigs))
		done <- true
	}()

	select {
	case <-done:
	case <-time.After(5 * time.Second):
		t.Fatal("Assembler deadlocked (MaxStepsPerWalk not working!)")
	}
}

func TestAssemblerMaxContigLengthEnforced(t *testing.T) {
	g := debruijn.NewGraph(4)
	longSeq := "ACGT"
	for i := 0; i < 500; i++ {
		longSeq += "ACGT"
	}
	_ = g.AddSequence(longSeq)

	asm := New()
	asm.MaxContigLength = 100
	contigs := asm.Assemble(g, 10)
	for _, c := range contigs {
		if c.Length > asm.MaxContigLength+10 {
			t.Errorf("Contig length %d exceeds max %d by more than k-1",
				c.Length, asm.MaxContigLength)
		}
	}
}

func TestAssemblerSnapshotDoesNotAffectOriginalGraph(t *testing.T) {
	g := debruijn.NewGraph(4)
	_ = g.AddSequence("ACGTGACTAGCT")

	nodesBefore := g.NodeCount()
	edgesBefore := g.EdgeCount()

	_ = Assemble(g, 5)

	if g.NodeCount() != nodesBefore {
		t.Errorf("Node count changed after assembly: %d -> %d",
			nodesBefore, g.NodeCount())
	}
	if g.EdgeCount() != edgesBefore {
		t.Errorf("Edge count changed after assembly: %d -> %d",
			edgesBefore, g.EdgeCount())
	}
}
