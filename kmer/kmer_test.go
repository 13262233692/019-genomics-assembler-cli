package kmer

import (
	"testing"
)

func TestExtract(t *testing.T) {
	seq := "ABCDEFGHIJ"
	kmers := Extract(seq, 3)
	if len(kmers) != 0 {
		t.Errorf("Expected 0 kmers for invalid bases, got %d", len(kmers))
	}

	seq = "ACGTACGT"
	kmers = Extract(seq, 3)
	expected := []string{"ACG", "CGT", "GTA", "TAC", "ACG", "CGT"}
	if len(kmers) != len(expected) {
		t.Fatalf("Expected %d kmers, got %d", len(expected), len(kmers))
	}
	for i, k := range kmers {
		if k != expected[i] {
			t.Errorf("Kmer %d: expected %s, got %s", i, expected[i], k)
		}
	}
}

func TestIsValid(t *testing.T) {
	tests := []struct {
		input    string
		expected bool
	}{
		{"ACGT", true},
		{"acgt", true},
		{"ACGTN", false},
		{"", true},
		{"A", true},
		{"X", false},
	}
	for _, tt := range tests {
		if IsValid(tt.input) != tt.expected {
			t.Errorf("IsValid(%s) = %v, expected %v", tt.input, !tt.expected, tt.expected)
		}
	}
}

func TestHash(t *testing.T) {
	h1 := Hash("ACGT")
	h2 := Hash("ACGT")
	h3 := Hash("TGCA")
	if h1 != h2 {
		t.Error("Same kmer should have same hash")
	}
	if h1 == h3 {
		t.Error("Different kmers should have different hashes")
	}
}

func TestReverseComplement(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"ACGT", "ACGT"},
		{"A", "T"},
		{"GCGC", "GCGC"},
		{"ATCG", "CGAT"},
		{"acgt", "acgt"},
		{"ATCGAT", "ATCGAT"},
	}
	for _, tt := range tests {
		result := ReverseComplement(tt.input)
		if result != tt.expected {
			t.Errorf("ReverseComplement(%s) = %s, expected %s", tt.input, result, tt.expected)
		}
	}
}

func TestCanonical(t *testing.T) {
	c1 := Canonical("ATCG")
	c2 := Canonical("CGAT")
	if c1 != c2 {
		t.Errorf("Canonical forms should be equal: %s vs %s", c1, c2)
	}
	if c1 != "ATCG" {
		t.Errorf("Expected canonical form 'ATCG', got '%s'", c1)
	}
}
