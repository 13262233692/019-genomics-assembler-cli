package debruijn

import (
	"sync"

	"github.com/genomics-assembler/cli/kmer"
)

type Node struct {
	Label    string
	Hash     uint64
	OutEdges []string
	InEdges  []string
}

type Graph struct {
	nodes map[string]*Node
	mu    sync.RWMutex
	k     int
}

func NewGraph(k int) *Graph {
	return &Graph{
		nodes: make(map[string]*Node),
		k:     k,
	}
}

func (g *Graph) K() int {
	return g.k
}

func (g *Graph) AddKmer(kmerStr string) {
	if len(kmerStr) != g.k {
		return
	}
	prefix := kmerStr[:g.k-1]
	suffix := kmerStr[1:]

	g.mu.Lock()
	defer g.mu.Unlock()

	if _, exists := g.nodes[prefix]; !exists {
		g.nodes[prefix] = &Node{
			Label: prefix,
			Hash:  kmer.Hash(prefix),
		}
	}
	if _, exists := g.nodes[suffix]; !exists {
		g.nodes[suffix] = &Node{
			Label: suffix,
			Hash:  kmer.Hash(suffix),
		}
	}

	prefixNode := g.nodes[prefix]
	suffixNode := g.nodes[suffix]

	prefixNode.OutEdges = appendUnique(prefixNode.OutEdges, suffix)
	suffixNode.InEdges = appendUnique(suffixNode.InEdges, prefix)
}

func (g *Graph) AddSequence(seq string) {
	kmers := kmer.Extract(seq, g.k)
	for _, k := range kmers {
		g.AddKmer(k)
	}
}

func (g *Graph) NodeCount() int {
	g.mu.RLock()
	defer g.mu.RUnlock()
	return len(g.nodes)
}

func (g *Graph) GetNode(key string) (*Node, bool) {
	g.mu.RLock()
	defer g.mu.RUnlock()
	n, ok := g.nodes[key]
	return n, ok
}

func (g *Graph) AllNodes() map[string]*Node {
	g.mu.RLock()
	defer g.mu.RUnlock()
	cp := make(map[string]*Node, len(g.nodes))
	for k, v := range g.nodes {
		cp[k] = v
	}
	return cp
}

func (g *Graph) OutDegree(key string) int {
	g.mu.RLock()
	defer g.mu.RUnlock()
	if n, ok := g.nodes[key]; ok {
		return len(n.OutEdges)
	}
	return 0
}

func (g *Graph) InDegree(key string) int {
	g.mu.RLock()
	defer g.mu.RUnlock()
	if n, ok := g.nodes[key]; ok {
		return len(n.InEdges)
	}
	return 0
}

func appendUnique(slice []string, val string) []string {
	for _, v := range slice {
		if v == val {
			return slice
		}
	}
	return append(slice, val)
}
