package restriction

import (
	"strings"
	"testing"
)

func TestFindEnzyme(t *testing.T) {
	tests := []struct {
		name     string
		expected string
		found    bool
	}{
		{"DLE-1", "DLE-1", true},
		{"BspQI", "BspQI", true},
		{"BbvCI", "BbvCI", true},
		{"unknown", "", false},
	}

	for _, tt := range tests {
		enzyme, ok := FindEnzyme(tt.name)
		if ok != tt.found {
			t.Errorf("FindEnzyme(%s): expected found=%v, got %v", tt.name, tt.found, ok)
		}
		if ok && enzyme.Name != tt.expected {
			t.Errorf("FindEnzyme(%s): expected name %s, got %s", tt.name, tt.expected, enzyme.Name)
		}
	}
}

func TestDigestSequence(t *testing.T) {
	seq := "NNNCTTAAGNNNNNNCTTAAGNNNNN"
	d := DigestSequence(1, seq, DLE1)

	if d == nil {
		t.Fatal("DigestSequence returned nil")
	}

	if d.ContigID != 1 {
		t.Errorf("Expected ContigID 1, got %d", d.ContigID)
	}

	if d.Length != len(seq) {
		t.Errorf("Expected length %d, got %d", len(seq), d.Length)
	}

	if d.SiteCount() < 2 {
		t.Errorf("Expected at least 2 sites, got %d", d.SiteCount())
	}

	if d.FirstSitePos() < 3 {
		t.Errorf("Expected first site >= 3, got %f", d.FirstSitePos())
	}

	if d.LastSitePos() <= d.FirstSitePos() {
		t.Errorf("Last site should be after first site")
	}
}

func TestDigestSequenceOrientation(t *testing.T) {
	seq := "NNNCTTAAGNNNN"
	d := DigestSequence(1, seq, DLE1)

	if d.SiteCount() != 1 {
		t.Errorf("Expected 1 site, got %d", d.SiteCount())
	}

	if d.Sites[0].Position < 3 || d.Sites[0].Position > 4 {
		t.Errorf("Expected position near 3, got %d", d.Sites[0].Position)
	}
}

func TestCompareDigests(t *testing.T) {
	seq1 := "NNNCTTAAGNNNCTTAAGNNNCTTAAG"
	seq2 := "NNNCTTAAGNNNCTTAAGNNNCTTAAG"

	d1 := DigestSequence(1, seq1, DLE1)
	d2 := DigestSequence(2, seq2, DLE1)

	similarity := CompareDigests(d1, d2, 10.0)
	if similarity <= 0.5 {
		t.Errorf("Expected high similarity for identical sequences, got %f", similarity)
	}

	seq3 := "NNNCTTAAGNNNNNNNNNNNNNNNNNNCTTAAG"
	d3 := DigestSequence(3, seq3, DLE1)
	similarity = CompareDigests(d1, d3, 10.0)
	if similarity >= 0.9 {
		t.Errorf("Expected lower similarity for different sequences, got %f", similarity)
	}
}

func TestEdgeKmers(t *testing.T) {
	longSeq := strings.Repeat("A", 100) + "CTTAAG" + strings.Repeat("G", 100)
	d := DigestSequence(1, longSeq, DLE1)

	if len(d.EdgeKmersL) == 0 {
		t.Error("Expected EdgeKmersL to be non-empty")
	}

	if len(d.EdgeKmersR) == 0 {
		t.Error("Expected EdgeKmersR to be non-empty")
	}
}

func TestDistances(t *testing.T) {
	seq := "NNNCTTAAGNNNCTTAAGNNN"
	d := DigestSequence(1, seq, DLE1)

	if len(d.Distances) != len(d.Positions)+1 {
		t.Errorf("Expected %d distances, got %d", len(d.Positions)+1, len(d.Distances))
	}

	total := 0.0
	for _, dist := range d.Distances {
		total += dist
	}
	if total != float64(d.Length) {
		t.Errorf("Sum of distances should equal length: %f vs %d", total, d.Length)
	}
}
