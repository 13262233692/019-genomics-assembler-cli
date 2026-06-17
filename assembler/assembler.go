package assembler

import (
	"strings"

	"github.com/genomics-assembler/cli/debruijn"
)

type Contig struct {
	ID       int
	Sequence string
	Length   int
}

func Assemble(g *debruijn.Graph, minLength int) []Contig {
	nodes := g.AllNodes()
	visited := make(map[string]bool)
	var contigs []Contig
	contigID := 0

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
			seq := walkForward(key, nodes, visited)
			if len(seq) >= minLength {
				contigID++
				contigs = append(contigs, Contig{
					ID:       contigID,
					Sequence: seq,
					Length:   len(seq),
				})
			}
		}
	}

	for key := range nodes {
		if visited[key] {
			continue
		}
		seq := walkForward(key, nodes, visited)
		if len(seq) >= minLength {
			contigID++
			contigs = append(contigs, Contig{
				ID:       contigID,
				Sequence: seq,
				Length:   len(seq),
			})
		}
	}

	return contigs
}

func walkForward(start string, nodes map[string]*debruijn.Node, visited map[string]bool) string {
	var builder strings.Builder
	current := start
	builder.WriteString(current)
	visited[current] = true

	for {
		node, ok := nodes[current]
		if !ok {
			break
		}
		if len(node.OutEdges) == 0 {
			break
		}
		next := node.OutEdges[0]
		if visited[next] {
			break
		}
		builder.WriteByte(next[len(next)-1])
		visited[next] = true
		current = next
	}

	return builder.String()
}
