package scaffold

import (
	"strings"
	"testing"

	"github.com/genomics-assembler/cli/cmap"
	"github.com/genomics-assembler/cli/restriction"
)

const testCMAP = `# CMAP File Version: 1.0
# EnzymeName: DLE-1
#
# CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
1	100000	5	1	1	10000	10.5	20	1
1	100000	5	2	1	30000	8.2	25	1
1	100000	5	3	1	55000	12.1	22	1
1	100000	5	4	1	75000	9.8	28	1
1	100000	5	5	1	95000	11.3	24	1
`

func TestDPAlignSites(t *testing.T) {
	query := []float64{10000, 30000, 55000}
	target := []float64{10000, 30000, 55000, 75000, 95000}

	score, qPath, tPath := DPAlignSites(query, target, 500.0)

	if score <= 0 {
		t.Errorf("Expected positive score, got %f", score)
	}

	if len(qPath) < 3 {
		t.Errorf("Expected at least 3 matches, got %d", len(qPath))
	}

	if len(qPath) != len(tPath) {
		t.Errorf("Path lengths should match: %d vs %d", len(qPath), len(tPath))
	}
}

func TestDPAlignSitesTolerance(t *testing.T) {
	query := []float64{10000, 30000}
	target := []float64{10050, 29950}

	score, _, _ := DPAlignSites(query, target, 100.0)
	if score <= 0 {
		t.Errorf("Expected positive score within tolerance, got %f", score)
	}

	score2, _, _ := DPAlignSites(query, target, 10.0)
	if score2 > score {
		t.Errorf("Tighter tolerance should give lower or equal score: %f vs %f", score2, score)
	}
}

func TestDPAlignSitesReverse(t *testing.T) {
	query := []float64{10000, 30000, 50000}
	target := []float64{0, 20000, 40000}

	revQuery := reverseFloat64(query, 50000)
	score, _, _ := DPAlignSites(revQuery, target, 500.0)

	if score <= 0 {
		t.Errorf("Reverse complement should still align, got %f", score)
	}
}

func TestAlignContigToMap(t *testing.T) {
	cmapData, err := cmap.ParseCMAP(strings.NewReader(testCMAP))
	if err != nil {
		t.Fatalf("Failed to parse CMAP: %v", err)
	}

	cm := cmapData.ContigMaps["1"]

	seq := strings.Repeat("N", 10000) + "CTTAAG" +
		strings.Repeat("N", 20000) + "CTTAAG" +
		strings.Repeat("N", 25000) + "CTTAAG" +
		strings.Repeat("N", 20000) + "CTTAAG" +
		strings.Repeat("N", 20000) + "CTTAAG"

	d := restriction.DigestSequence(1, seq, restriction.DLE1)

	s := NewScaffolder()
	s.MinMatches = 3
	s.Tolerance = 500.0

	ali := s.AlignContigToMap(d, cm)
	if ali == nil {
		t.Fatal("Expected alignment, got nil")
	}

	if ali.Score <= 0 {
		t.Errorf("Expected positive score, got %f", ali.Score)
	}

	if ali.ContigID != 1 {
		t.Errorf("Expected ContigID 1, got %d", ali.ContigID)
	}
}

func TestBuildScaffold(t *testing.T) {
	cmapData, _ := cmap.ParseCMAP(strings.NewReader(testCMAP))
	cm := cmapData.ContigMaps["1"]

	seq1 := strings.Repeat("N", 10000) + "CTTAAG" + strings.Repeat("N", 20000)
	seq2 := strings.Repeat("N", 5000) + "CTTAAG" + strings.Repeat("N", 25000)

	d1 := restriction.DigestSequence(1, seq1, restriction.DLE1)
	d2 := restriction.DigestSequence(2, seq2, restriction.DLE1)

	digests := map[int]*restriction.ContigDigest{1: d1, 2: d2}

	ali1 := &Alignment{
		ContigID:  1,
		MapID:     "1",
		Score:     0.9,
		Orientation: 1,
		MapStartPos: 0,
		MapEndPos:   30000,
	}
	ali2 := &Alignment{
		ContigID:  2,
		MapID:     "1",
		Score:     0.9,
		Orientation: 1,
		MapStartPos: 30000,
		MapEndPos:   60000,
	}

	s := NewScaffolder()
	sc := s.BuildScaffold(1, []*Alignment{ali1, ali2}, digests, cm)

	if sc == nil {
		t.Fatal("BuildScaffold returned nil")
	}

	if sc.ContigCount != 2 {
		t.Errorf("Expected 2 contigs, got %d", sc.ContigCount)
	}

	if sc.TotalGapBases <= 0 {
		t.Error("Expected gap bases between contigs")
	}

	if sc.Length != len(seq1)+len(seq2)+sc.TotalGapBases {
		t.Errorf("Scaffold length mismatch: expected %d, got %d",
			len(seq1)+len(seq2)+sc.TotalGapBases, sc.Length)
	}

	if !strings.Contains(sc.Sequence, "N") {
		t.Error("Expected N gaps in scaffold sequence")
	}
}

func TestScaffolderRun(t *testing.T) {
	cmapData, _ := cmap.ParseCMAP(strings.NewReader(testCMAP))

	seq1 := strings.Repeat("N", 10000) + "CTTAAG" +
		strings.Repeat("N", 20000) + "CTTAAG" +
		strings.Repeat("N", 25000)
	seq2 := strings.Repeat("N", 1000) + "CTTAAG" +
		strings.Repeat("N", 19000) + "CTTAAG" +
		strings.Repeat("N", 20000)

	d1 := restriction.DigestSequence(1, seq1, restriction.DLE1)
	d2 := restriction.DigestSequence(2, seq2, restriction.DLE1)
	d3 := restriction.DigestSequence(3, strings.Repeat("A", 1000), restriction.DLE1)

	digests := []*restriction.ContigDigest{d1, d2, d3}

	s := NewScaffolder()
	s.MinMatches = 2
	s.Tolerance = 2000.0
	s.MinScore = 0.1

	scaffolds := s.Run(digests, cmapData)

	if len(scaffolds) == 0 {
		t.Fatal("Expected at least one scaffold")
	}

	mappedCount := 0
	for _, sc := range scaffolds {
		if !strings.HasPrefix(sc.MapID, "unmapped") {
			mappedCount++
		}
	}

	if mappedCount == 0 {
		t.Error("Expected at least one mapped scaffold")
	}
}

func TestScaffoldOrientation(t *testing.T) {
	seq := "ATGCATGCCTTAAGCGTACGTA"
	d := restriction.DigestSequence(1, seq, restriction.DLE1)

	revSeq := reverseComplementSeq(seq)
	dRev := restriction.DigestSequence(2, revSeq, restriction.DLE1)

	aliFwd := &Alignment{
		ContigID:    1,
		Orientation: 1,
		MapStartPos: 0,
		MapEndPos:   float64(len(seq)),
	}
	aliRev := &Alignment{
		ContigID:    2,
		Orientation: -1,
		MapStartPos: float64(len(seq)),
		MapEndPos:   float64(len(seq) * 2),
	}

	digests := map[int]*restriction.ContigDigest{1: d, 2: dRev}

	cmapData, _ := cmap.ParseCMAP(strings.NewReader(testCMAP))
	cm := cmapData.ContigMaps["1"]

	s := NewScaffolder()
	sc := s.BuildScaffold(1, []*Alignment{aliFwd, aliRev}, digests, cm)

	if sc == nil {
		t.Fatal("BuildScaffold returned nil")
	}

	if len(sc.Sequence) <= len(seq) {
		t.Error("Scaffold should contain both contigs and gap")
	}
}

func TestScaffolderDefaults(t *testing.T) {
	s := NewScaffolder()

	if s.Tolerance != DefaultTolerance {
		t.Errorf("Default tolerance mismatch")
	}
	if s.MinMatches != DefaultMinMatches {
		t.Errorf("Default min matches mismatch")
	}
	if s.MinScore != DefaultMinScore {
		t.Errorf("Default min score mismatch")
	}
	if s.GapSize != DefaultGapSize {
		t.Errorf("Default gap size mismatch")
	}
}

func TestReverseFloat64(t *testing.T) {
	input := []float64{10000, 30000, 50000}
	total := 50000.0
	result := reverseFloat64(input, total)

	expected := []float64{0, 20000, 40000}
	for i, v := range result {
		if v != expected[i] {
			t.Errorf("reverseFloat64[%d]: expected %f, got %f", i, expected[i], v)
		}
	}
}
