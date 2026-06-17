package scaffold

import (
	"fmt"
	"log"
	"math"
	"sort"
	"strings"

	"github.com/genomics-assembler/cli/cmap"
	"github.com/genomics-assembler/cli/restriction"
)

const (
	DefaultTolerance    = 500.0
	DefaultMinMatches   = 3
	DefaultMinScore  = 0.3
	DefaultGapSize = 500
)

type Alignment struct {
	ContigID    int
	MapID       string
	Score       float64
	StartIdx    int
	EndIdx      int
	MapStartIdx int
	MapEndIdx   int
	Orientation int
	StartPos    float64
	EndPos      float64
	MapStartPos float64
	MapEndPos   float64
}

type Scaffold struct {
	ID          int
	MapID        string
	Length       int
	Sequence     string
	Components   []ScaffoldComponent
	TotalGapBases int
	ContigCount  int
}

type ScaffoldComponent struct {
	ContigID    int
	Orientation int
	StartPos    int
	EndPos      int
	GapAfter    int
}

type Scaffolder struct {
	Tolerance   float64
	MinMatches  int
	MinScore    float64
	GapSize     int
	Verbose     bool
}

func NewScaffolder() *Scaffolder {
	return &Scaffolder{
		Tolerance:  DefaultTolerance,
		MinMatches: DefaultMinMatches,
		MinScore:   DefaultMinScore,
		GapSize:    DefaultGapSize,
	}
}

func DPAlignSites(query, target []float64, tolerance float64) (float64, []int, []int) {
	m := len(query)
	n := len(target)
	if m == 0 || n == 0 {
		return 0, nil, nil
	}

	dp := make([][]float64, m+1)
	for i := range dp {
		dp[i] = make([]float64, n+1)
	}

	backtrack := make([][]int, m+1)
	for i := range backtrack {
		backtrack[i] = make([]int, n+1)
	}

	gap := -1.0

	for i := 1; i <= m; i++ {
		for j := 1; j <= n; j++ {
			diff := math.Abs(query[i-1] - target[j-1])
			var match float64
			if diff <= tolerance {
				match = 1.0 / (1.0 + diff/tolerance)
			} else {
				match = -2.0
			}

			diag := dp[i-1][j-1] + match
			up := dp[i-1][j] + gap
			left := dp[i][j-1] + gap

			if diag >= up && diag >= left {
				dp[i][j] = diag
				backtrack[i][j] = 1
			} else if up >= left {
				dp[i][j] = up
				backtrack[i][j] = 2
			} else {
				dp[i][j] = left
				backtrack[i][j] = 3
			}
		}
	}

	bestScore := 0.0
	bestI, bestJ := 0, 0
	for i := 0; i <= m; i++ {
		for j := 0; j <= n; j++ {
			if dp[i][j] > bestScore {
				bestScore = dp[i][j]
				bestI = i
				bestJ = j
			}
		}
	}

	queryPath := make([]int, 0)
	targetPath := make([]int, 0)
	i, j := bestI, bestJ

	for i > 0 && j > 0 && dp[i][j] > 0 {
		switch backtrack[i][j] {
		case 1:
			i--
			j--
			queryPath = append(queryPath, i)
			targetPath = append(targetPath, j)
		case 2:
			i--
		case 3:
			j--
		default:
			break
		}
	}

	reverseInts(queryPath)
	reverseInts(targetPath)

	return bestScore, queryPath, targetPath
}

func (s *Scaffolder) AlignContigToMap(
	digest *restriction.ContigDigest, cm *cmap.ContigMap) *Alignment {
	if digest == nil || cm == nil {
		return nil
	}
	if digest.SiteCount() < s.MinMatches || len(cm.SitePositions) < s.MinMatches {
		return nil
	}

	scoreFwd, qPathFwd, tPathFwd := DPAlignSites(digest.Positions, cm.SitePositions, s.Tolerance)

	revPositions := reverseFloat64(digest.Positions, float64(digest.Length))
	scoreRev, qPathRev, tPathRev := DPAlignSites(revPositions, cm.SitePositions, s.Tolerance)

	score := scoreFwd
	orientation := 1
	qPath := qPathFwd
	tPath := tPathFwd
	startPos := digest.FirstSitePos()
	endPos := digest.LastSitePos()

	if scoreRev > scoreFwd {
		score = scoreRev
		orientation = -1
		qPath = qPathRev
		tPath = tPathRev
		startPos = float64(digest.Length) - digest.LastSitePos()
		endPos = float64(digest.Length) - digest.FirstSitePos()
	}

	maxPossible := float64(min(len(digest.Positions), len(cm.SitePositions)))
	normScore := score / maxPossible

	if len(qPath) < s.MinMatches || normScore < s.MinScore {
		return nil
	}

	mapStartIdx := 0
	mapEndIdx := len(cm.SitePositions) - 1
	if len(tPath) > 0 {
		mapStartIdx = tPath[0] - 1
		if mapStartIdx < 0 {
			mapStartIdx = 0
		}
		mapEndIdx = tPath[len(tPath)-1] - 1
		if mapEndIdx >= len(cm.SitePositions) {
			mapEndIdx = len(cm.SitePositions) - 1
		}
	}

	mapStartPos := 0.0
	mapEndPos := cm.Length
	if mapStartIdx < len(cm.SitePositions) {
		mapStartPos = cm.SitePositions[mapStartIdx]
	}
	if mapEndIdx >= 0 && mapEndIdx < len(cm.SitePositions) {
		mapEndPos = cm.SitePositions[mapEndIdx]
	}

	return &Alignment{
		ContigID:    digest.ContigID,
		MapID:       cm.ID,
		Score:         normScore,
		StartIdx:      qPath[0],
		EndIdx:        qPath[len(qPath)-1],
		MapStartIdx:   mapStartIdx,
		MapEndIdx:     mapEndIdx,
		Orientation:   orientation,
		StartPos:      startPos,
		EndPos:        endPos,
		MapStartPos:   mapStartPos,
		MapEndPos:     mapEndPos,
	}
}

func (s *Scaffolder) OrderContigsOnMap(
	digests []*restriction.ContigDigest,
	cm *cmap.ContigMap) []*Alignment {

	alignments := make([]*Alignment, 0, len(digests))

	for _, d := range digests {
		ali := s.AlignContigToMap(d, cm)
		if ali != nil {
			alignments = append(alignments, ali)
		}
	}

	usedContigs := make(map[int]bool)
	final := make([]*Alignment, 0)

	sort.Slice(alignments, func(i, j int) bool {
		return alignments[i].MapStartPos < alignments[j].MapStartPos
	})

	for _, ali := range alignments {
		if usedContigs[ali.ContigID] {
			continue
		}
		conflict := false
		for _, f := range final {
			overlap := !(ali.MapEndPos < f.MapStartPos-1000 || ali.MapStartPos > f.MapEndPos+1000)
			if overlap && ali.ContigID != f.ContigID {
				conflict = true
				break
			}
		}
		if !conflict {
			final = append(final, ali)
			usedContigs[ali.ContigID] = true
		}
	}

	sort.Slice(final, func(i, j int) bool {
		return final[i].MapStartPos < final[j].MapStartPos
	})

	return final
}

func (s *Scaffolder) BuildScaffold(
	scaffoldID int,
	alignments []*Alignment,
	digests map[int]*restriction.ContigDigest,
	cm *cmap.ContigMap) *Scaffold {

	if len(alignments) == 0 {
		return nil
	}

	sc := &Scaffold{
		ID:         scaffoldID,
		MapID:      cm.ID,
		Components: make([]ScaffoldComponent, 0, len(alignments)),
	}

	var builder strings.Builder
	currentPos := 0

	for i, ali := range alignments {
		d := digests[ali.ContigID]
		if d == nil {
			continue
		}

		seq := d.Sequence
		if ali.Orientation == -1 {
			seq = reverseComplementSeq(seq)
		}

		gap := 0
		if i > 0 {
			prevAli := alignments[i-1]
			prevD := digests[prevAli.ContigID]
			if prevD != nil {
				estimatedGap := int(ali.MapStartPos - prevAli.MapEndPos)
				if estimatedGap > 0 {
					gap = estimatedGap
				} else {
					gap = s.GapSize
				}
			} else {
				gap = s.GapSize
			}
			if gap > 100000 {
				gap = 100000
			}
			if gap < 100 {
				gap = s.GapSize
			}
			builder.WriteString(strings.Repeat("N", gap))
			sc.TotalGapBases += gap
			currentPos += gap
		}

		start := currentPos
		builder.WriteString(seq)
		end := currentPos + len(seq)

		sc.Components = append(sc.Components, ScaffoldComponent{
			ContigID:    ali.ContigID,
			Orientation: ali.Orientation,
			StartPos:    start,
			EndPos:      end,
			GapAfter:    gap,
		})

		currentPos = end
		sc.ContigCount++
	}

	sc.Sequence = builder.String()
	sc.Length = len(sc.Sequence)

	return sc
}

func (s *Scaffolder) Run(
	digests []*restriction.ContigDigest,
	cmapData *cmap.CMAP) []*Scaffold {

	if s.Verbose {
		log.Printf("Scaffolder: %d contigs to align to %d maps",
			len(digests), cmapData.ContigCount())
	}

	digestMap := make(map[int]*restriction.ContigDigest)
	for _, d := range digests {
		digestMap[d.ContigID] = d
	}

	scaffolds := make([]*Scaffold, 0)
	scaffoldID := 1

	for _, mapID := range cmapData.OrderedIDs {
		cm := cmapData.ContigMaps[mapID]
		if cm == nil || len(cm.SitePositions) < s.MinMatches {
			continue
		}

		alignments := s.OrderContigsOnMap(digests, cm)
		if len(alignments) < 1 {
			continue
		}

		sc := s.BuildScaffold(scaffoldID, alignments, digestMap, cm)
		if sc != nil && sc.Length > 0 {
			scaffolds = append(scaffolds, sc)
			scaffoldID++
			if s.Verbose {
				log.Printf("Scaffold %d (map %s): %d contigs, length %d",
					sc.ID, sc.MapID, sc.ContigCount, sc.Length)
			}
		}
	}

	unmappedDigests := make([]*restriction.ContigDigest, 0)
	for _, d := range digests {
		mapped := false
		for _, sc := range scaffolds {
			for _, comp := range sc.Components {
				if comp.ContigID == d.ContigID {
					mapped = true
					break
				}
			}
		}
		if !mapped {
			unmappedDigests = append(unmappedDigests, d)
		}
	}

	for _, d := range unmappedDigests {
		sc := &Scaffold{
			ID:         scaffoldID,
			MapID:      fmt.Sprintf("unmapped_%d", scaffoldID),
			Length:     d.Length,
			Sequence:   d.Sequence,
			Components: []ScaffoldComponent{{
				ContigID:    d.ContigID,
				Orientation: 1,
				StartPos:    0,
				EndPos:      d.Length,
			}},
			ContigCount: 1,
		}
		scaffolds = append(scaffolds, sc)
		scaffoldID++
	}

	if s.Verbose {
		log.Printf("Scaffolder complete: %d scaffolds built (%d mapped, %d unmapped)",
			len(scaffolds), len(scaffolds)-len(unmappedDigests), len(unmappedDigests))
	}

	return scaffolds
}

func reverseInts(a []int) {
	for i, j := 0, len(a)-1; i < j; i, j = i+1, j-1 {
		a[i], a[j] = a[j], a[i]
	}
}

func reverseFloat64(a []float64, totalLen float64) []float64 {
	result := make([]float64, len(a))
	for i, v := range a {
		result[len(a)-1-i] = totalLen - v
	}
	return result
}

func reverseComplementSeq(s string) string {
	comp := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
		'N': 'N', 'n': 'n',
		'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
	}
	result := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		result[len(s)-1-i] = comp[s[i]]
	}
	return string(result)
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
