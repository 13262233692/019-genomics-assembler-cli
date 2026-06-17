package assembler

import (
	"log"
	"strings"

	"github.com/genomics-assembler/cli/debruijn"
)

const (
	DefaultMaxContigLength = 10_000_000
	DefaultMaxStepsPerWalk     = 500_000_000
	ProbeInterval        = 1_000_000
)

type Contig struct {
	ID       int
	Sequence string
	Length   int
}

type Assembler struct {
	MaxContigLength int
	MaxStepsPerWalk int
	Verbose         bool
}

func New() *Assembler {
	return &Assembler{
		MaxContigLength: DefaultMaxContigLength,
		MaxStepsPerWalk: DefaultMaxStepsPerWalk,
	}
}

func Assemble(g *debruijn.Graph, minLength int) []Contig {
	a := New()
	return a.Assemble(g, minLength)
}

func (a *Assembler) Assemble(g *debruijn.Graph, minLength int) []Contig {
	snap := g.Snapshot()
	nodes := snap.Nodes
	totalNodes := len(nodes)

	if a.Verbose {
		log.Printf("Assembler: snapshot with %d nodes, starting traversal", totalNodes)
	}

	visited := make(map[string]bool, totalNodes)
	var contigs []Contig
	contigID := 0
	stepsUsed := uint64(0)
	deadlocks := 0
	cycles := 0

	for key, node := range nodes {
		if visited[key] {
			continue
		}
		if len(node.InEdges) == 0 && len(node.OutEdges) == 0 {
			visited[key] = true
			seq := key
			if len(seq) >= minLength {
				contigID++
				contigs = append(contigs, Contig{
					ID:       contigID,
					Sequence: seq,
					Length:   len(seq),
				})
			}
			continue
		}
		if len(node.InEdges) == 0 {
			result := a.walkSafe(key, nodes, visited, &stepsUsed, &deadlocks, &cycles)
			if len(result) >= minLength {
				contigID++
				contigs = append(contigs, Contig{
					ID:       contigID,
					Sequence: result,
					Length:   len(result),
				})
			}
		}
	}

	for key := range nodes {
		if visited[key] {
			continue
		}
		result := a.walkSafe(key, nodes, visited, &stepsUsed, &deadlocks, &cycles)
		if len(result) >= minLength {
			contigID++
			contigs = append(contigs, Contig{
				ID:       contigID,
				Sequence: result,
				Length:   len(result),
			})
		}
	}

	if a.Verbose {
		log.Printf("Assembler: total steps=%d, deadlock-breaks=%d, cycle-breaks=%d, contigs=%d",
			stepsUsed, deadlocks, cycles, contigID)
	}

	return contigs
}

func (a *Assembler) walkSafe(
	start string,
	nodes map[string]*debruijn.NodeSnapshot,
	visited map[string]bool,
	totalSteps *uint64,
	deadlockCount *int,
	cycleCount *int,
) string {
	var builder strings.Builder
	current := start
	builder.WriteString(current)
	visited[current] = true

	walkerSteps := 0
	tortoise := start
	tortoiseStep := 0
	cycleDetectEvery := 256

	for {
		if *totalSteps > uint64(a.MaxStepsPerWalk) {
			*deadlockCount++
			if a.Verbose {
				log.Printf("DEADLOCK PROBE: global step limit %d exceeded at node %q, aborting walk",
					a.MaxStepsPerWalk, current)
			}
			break
		}

		node, ok := nodes[current]
		if !ok {
			break
		}
		if len(node.OutEdges) == 0 {
			break
		}

		next := pickNextEdge(node.OutEdges, visited)
		if next == "" {
			break
		}

		if builder.Len() >= a.MaxContigLength {
			if a.Verbose {
				log.Printf("LENGTH PROBE: contig length %d hit cap, stopping at node %q",
					builder.Len(), current)
			}
			break
		}

		if visited[next] {
			break
		}

		builder.WriteByte(next[len(next)-1])
		visited[next] = true
		current = next
		*totalSteps++
		walkerSteps++

		if walkerSteps%cycleDetectEvery == 0 {
			if walkerSteps-tortoiseStep >= cycleDetectEvery {
				if tortoise == current {
					*cycleCount++
					if a.Verbose {
						log.Printf("CYCLE PROBE: Floyd tortoise-hare cycle detected starting at %q (step=%d), breaking",
							current, walkerSteps)
					}
					break
				}
				tortoise = pickTortoiseNext(tortoise, nodes, visited)
				tortoiseStep = walkerSteps
			}
		}

		if walkerSteps%ProbeInterval == 0 && a.Verbose {
			log.Printf("PROBE: walker step=%d, contig_len=%d, node=%q",
				walkerSteps, builder.Len(), current)
		}
	}

	return builder.String()
}

func pickNextEdge(edges []string, visited map[string]bool) string {
	for _, e := range edges {
		if !visited[e] {
			return e
		}
	}
	for _, e := range edges {
		return e
	}
	return ""
}

func pickTortoiseNext(nodeKey string, nodes map[string]*debruijn.NodeSnapshot, visited map[string]bool) string {
	node, ok := nodes[nodeKey]
	if !ok {
		return nodeKey
	}
	if len(node.OutEdges) == 0 {
		return nodeKey
	}
	return node.OutEdges[0]
}
