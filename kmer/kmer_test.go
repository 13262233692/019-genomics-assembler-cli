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

func TestHash2(t *testing.T) {
	h1 := Hash2("ACGT")
	h2 := Hash2("ACGT")
	h3 := Hash2("TGCA")
	if h1 != h2 {
		t.Error("Same kmer should have same Hash2")
	}
	if h1 == h3 {
		t.Error("Different kmers should have different Hash2")
	}
}

func TestHashDistinctFromHash2(t *testing.T) {
	h1 := Hash("ACGTACGTAGCTAGCT")
	h2 := Hash2("ACGTACGTAGCTAGCT")
	if h1 == h2 {
		t.Log("Warning: Hash and Hash2 produced same value (unlikely but possible)")
	}
}

func TestDoubleHash(t *testing.T) {
	dh1 := HashDouble("AAAAAAAAAAAAAAA")
	dh2 := HashDouble("AAAAAAAAAAAAAAA")
	dh3 := HashDouble("TTTTTTTTTTTTTTT")
	if !dh1.Equals(dh2) {
		t.Error("Same kmer should have equal DoubleHash")
	}
	if dh1.Equals(dh3) {
		t.Error("Different kmers should have different DoubleHash")
	}
	if dh1.H1 == dh3.H1 && dh1.H2 == dh3.H2 {
		t.Error("Different kmers produced double-hash collision")
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

func TestEncodeDecodeKmer(t *testing.T) {
	tests := []string{"ACGT", "AAAAAAAA", "GCGCGCGC", "ATATATAT", "ACGTACGT"}
	for _, seq := range tests {
		if len(seq) > 32 {
			continue
		}
		encoded := EncodeKmer(seq)
		decoded := DecodeKmer(encoded, len(seq))
		if decoded != seq {
			t.Errorf("Encode/Decode roundtrip failed: %s -> %d -> %s", seq, encoded, decoded)
		}
	}
}

func TestEncodeKmerDistinct(t *testing.T) {
	seen := make(map[uint64]string)
	kmers := []string{"AAAA", "AAAC", "AAAG", "AAAT", "CAAA", "GAAA", "TAAA", "ACGT"}
	for _, k := range kmers {
		v := EncodeKmer(k)
		if prev, exists := seen[v]; exists {
			t.Errorf("Collision: %s and %s both encode to %d", prev, k, v)
		}
		seen[v] = k
	}
}
