package restriction

import (
	"strings"
)

type Enzyme struct {
	Name        string
	Recognition string
	Length      int
	CutOffset   int
}

var (
	BspQI  = Enzyme{Name: "BspQI", Recognition: "GCTCTTC", Length: 7, CutOffset: 1}
	BbvCI  = Enzyme{Name: "BbvCI", Recognition: "CCTCAGC", Length: 7, CutOffset: 2}
	BssSI  = Enzyme{Name: "BssSI", Recognition: "CACGAG", Length: 6, CutOffset: 1}
	DLE1   = Enzyme{Name: "DLE-1", Recognition: "CTTAAG", Length: 6, CutOffset: 1}
	SwaI   = Enzyme{Name: "SwaI", Recognition: "ATTTAAAT", Length: 8, CutOffset: 3}
	PacI   = Enzyme{Name: "PacI", Recognition: "TTAATTAA", Length: 8, CutOffset: 3}
)

type Site struct {
	Position int
	Strand   int
	Sequence string
}

type ContigDigest struct {
	ContigID     int
	Sequence     string
	Length       int
	Enzyme       Enzyme
	Sites        []Site
	Positions    []float64
	Distances    []float64
	EdgeKmersL   []string
	EdgeKmersR   []string
}

var EnzymeMap = map[string]Enzyme{
	"BspQI":   BspQI,
	"BbvCI":   BbvCI,
	"BssSI":   BssSI,
	"DLE-1":   DLE1,
	"DLE1":    DLE1,
	"SwaI":    SwaI,
	"PacI":    PacI,
}

func FindEnzyme(name string) (Enzyme, bool) {
	if e, ok := EnzymeMap[name]; ok {
		return e, true
	}
	for eName, e := range EnzymeMap {
		if strings.EqualFold(eName, name) {
			return e, true
		}
	}
	return Enzyme{}, false
}

func DigestSequence(contigID int, sequence string, enzyme Enzyme) *ContigDigest {
	seq := strings.ToUpper(sequence)
	seqLen := len(seq)
	sites := make([]Site, 0)
	positions := make([]float64, 0)
	motif := strings.ToUpper(enzyme.Recognition)
	motifLen := len(motif)

	for i := 0; i <= seqLen-motifLen; i++ {
		if seq[i:i+motifLen] == motif {
			site := Site{
				Position: i + enzyme.CutOffset,
				Strand:   1,
				Sequence: motif,
			}
			sites = append(sites, site)
			positions = append(positions, float64(site.Position))
		}
		revMotif := reverseComplement(motif)
		if motif != revMotif && seq[i:i+motifLen] == revMotif {
			site := Site{
				Position: i + motifLen - enzyme.CutOffset,
				Strand:   -1,
				Sequence: revMotif,
			}
			sites = append(sites, site)
			positions = append(positions, float64(site.Position))
		}
	}

	sortSitesByPosition(sites, positions)

	distances := make([]float64, 0, len(positions)+1)
	prev := 0.0
	for _, p := range positions {
		distances = append(distances, p-prev)
		prev = p
	}
	if seqLen > 0 {
		distances = append(distances, float64(seqLen)-prev)
	}

	edgeK := 31
	if seqLen < edgeK*2 {
		edgeK = seqLen / 2
	}
	edgeKmersL := make([]string, 0)
	edgeKmersR := make([]string, 0)
	if seqLen >= edgeK {
		for i := 0; i <= edgeK-15; i++ {
			if i+15 <= seqLen {
				edgeKmersL = append(edgeKmersL, seq[i:i+15])
			}
		}
		for i := seqLen - 15; i >= seqLen-edgeK; i-- {
			if i >= 0 && i+15 <= seqLen {
				edgeKmersR = append(edgeKmersR, seq[i:i+15])
			}
		}
	}

	return &ContigDigest{
		ContigID:   contigID,
		Sequence:   seq,
		Length:     seqLen,
		Enzyme:     enzyme,
		Sites:      sites,
		Positions:  positions,
		Distances:  distances,
		EdgeKmersL: edgeKmersL,
		EdgeKmersR: edgeKmersR,
	}
}

func (d *ContigDigest) SiteCount() int {
	return len(d.Sites)
}

func (d *ContigDigest) FirstSitePos() float64 {
	if len(d.Positions) > 0 {
		return d.Positions[0]
	}
	return -1
}

func (d *ContigDigest) LastSitePos() float64 {
	if len(d.Positions) > 0 {
		return d.Positions[len(d.Positions)-1]
	}
	return -1
}

func reverseComplement(s string) string {
	comp := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
	}
	result := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		result[len(s)-1-i] = comp[s[i]]
	}
	return string(result)
}

func sortSitesByPosition(sites []Site, positions []float64) {
	for i := 1; i < len(sites); i++ {
		for j := i; j > 0 && sites[j].Position < sites[j-1].Position; j-- {
			sites[j], sites[j-1] = sites[j-1], sites[j]
			positions[j], positions[j-1] = positions[j-1], positions[j]
		}
	}
}

func CompareDigests(a, b *ContigDigest, tolerance float64) float64 {
	if a == nil || b == nil {
		return 0
	}
	if a.Enzyme.Name != b.Enzyme.Name {
		return 0
	}
	matches := 0
	i, j := 0, 0
	for i < len(a.Positions) && j < len(b.Positions) {
		diff := a.Positions[i] - b.Positions[j]
		if diff < -tolerance {
			i++
		} else if diff > tolerance {
			j++
		} else {
			matches++
			i++
			j++
		}
	}
	total := float64(len(a.Positions) + len(b.Positions))
	if total == 0 {
		return 0
	}
	return 2.0 * float64(matches) / total
}
