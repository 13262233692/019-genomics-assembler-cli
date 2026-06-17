package assembler

import (
	"testing"

	"github.com/genomics-assembler/cli/debruijn"
)

func TestAssembleSimple(t *testing.T) {
	g := debruijn.NewGraph(4)
	g.AddSequence("ACGTGACTAGCT")

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
	g.AddSequence("ACGT")

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
	g.AddSequence("AAGCTTGC")
	g.AddSequence("CCGGAATT")

	contigs := Assemble(g, 5)
	if len(contigs) < 2 {
		for i, c := range contigs {
			t.Logf("Contig %d: len=%d, seq=%s", i, c.Length, c.Sequence)
		}
		t.Errorf("Expected at least 2 contigs, got %d", len(contigs))
	}
}
