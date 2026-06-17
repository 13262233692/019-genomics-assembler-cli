package cmap

import (
	"strings"
	"testing"
)

const testCMAP = `# CMAP File Version: 1.0
# EnzymeName: DLE-1
# EnzymeSequence: CTTAAG
#
# CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
1	100000	3	1	1	10000	10.5	20	1
1	100000	3	2	1	30000	8.2	25	1
1	100000	3	3	1	50000	12.1	22	1
2	50000	2	1	1	15000	7.5	18	1
2	50000	2	2	1	40000	9.0	22	1
`

func TestParseCMAP(t *testing.T) {
	cmap, err := ParseCMAP(strings.NewReader(testCMAP))
	if err != nil {
		t.Fatalf("ParseCMAP failed: %v", err)
	}

	if cmap.EnzymeName != "DLE-1" {
		t.Errorf("Expected enzyme DLE-1, got %s", cmap.EnzymeName)
	}

	if cmap.ContigCount() != 2 {
		t.Errorf("Expected 2 contigs, got %d", cmap.ContigCount())
	}

	if cmap.TotalSites() != 5 {
		t.Errorf("Expected 5 total sites, got %d", cmap.TotalSites())
	}

	cm, ok := cmap.GetContig("1")
	if !ok {
		t.Fatal("Contig 1 not found")
	}

	if cm.Length != 100000 {
		t.Errorf("Expected length 100000, got %f", cm.Length)
	}

	if len(cm.Sites) != 3 {
		t.Errorf("Expected 3 sites, got %d", len(cm.Sites))
	}

	if cm.Sites[0].Position != 10000 {
		t.Errorf("Expected first site at 10000, got %f", cm.Sites[0].Position)
	}

	if cm.Sites[1].Position != 30000 {
		t.Errorf("Expected second site at 30000, got %f", cm.Sites[1].Position)
	}

	expectedPositions := []float64{10000, 30000, 50000}
	for i, p := range cm.SitePositions {
		if p != expectedPositions[i] {
			t.Errorf("SitePositions[%d]: expected %f, got %f", i, expectedPositions[i], p)
		}
	}
}

func TestLongestContig(t *testing.T) {
	cmap, err := ParseCMAP(strings.NewReader(testCMAP))
	if err != nil {
		t.Fatalf("ParseCMAP failed: %v", err)
	}

	longest := cmap.LongestContig()
	if longest == nil {
		t.Fatal("LongestContig returned nil")
	}
	if longest.ID != "1" {
		t.Errorf("Expected longest contig 1, got %s", longest.ID)
	}
}

func TestOpenCMAP(t *testing.T) {
	cmap, err := OpenCMAP("../testdata/test.cmap")
	if err != nil {
		t.Fatalf("OpenCMAP failed: %v", err)
	}

	if cmap.ContigCount() < 3 {
		t.Errorf("Expected at least 3 contigs, got %d", cmap.ContigCount())
	}

	if cmap.TotalSites() < 10 {
		t.Errorf("Expected at least 10 sites, got %d", cmap.TotalSites())
	}
}
