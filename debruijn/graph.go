package debruijn

import (
	"errors"
	"log"
	"sync"
	"sync/atomic"

	"github.com/genomics-assembler/cli/kmer"
)

const (
	DefaultMaxNodes   = 500_000_000
	DefaultMaxEdges   = 2_000_000_000
	EdgeBucketInitial = 2
)

var (
	ErrNodeLimitExceeded = errors.New("maximum node limit exceeded")
	ErrEdgeLimitExceeded = errors.New("maximum edge limit exceeded")
)

type Node struct {
	Label    string
	Hash     kmer.DoubleHash
	OutEdges map[string]struct{}
	InEdges  map[string]struct{}
	RefCount uint32
}

type Graph struct {
	nodes     map[string]*Node
	mu        sync.RWMutex
	k         int
	nodeCount uint64
	edgeCount uint64
	maxNodes  uint64
	maxEdges  uint64
}

func NewGraph(k int) *Graph {
	return &Graph{
		nodes:    make(map[string]*Node),
		k:        k,
		maxNodes: DefaultMaxNodes,
		maxEdges: DefaultMaxEdges,
	}
}

func NewGraphWithLimits(k int, maxNodes, maxEdges uint64) *Graph {
	return &Graph{
		nodes:    make(map[string]*Node),
		k:        k,
		maxNodes: maxNodes,
		maxEdges: maxEdges,
	}
}

func (g *Graph) K() int {
	return g.k
}

func (g *Graph) NodeCount() int {
	return int(atomic.LoadUint64(&g.nodeCount))
}

func (g *Graph) EdgeCount() uint64 {
	return atomic.LoadUint64(&g.edgeCount)
}

func (g *Graph) SetMaxNodes(n uint64) {
	atomic.StoreUint64(&g.maxNodes, n)
}

func (g *Graph) SetMaxEdges(n uint64) {
	atomic.StoreUint64(&g.maxEdges, n)
}

func (g *Graph) AddKmer(kmerStr string) error {
	if len(kmerStr) != g.k {
		return nil
	}
	prefix := kmerStr[:g.k-1]
	suffix := kmerStr[1:]

	g.mu.Lock()
	defer g.mu.Unlock()

	if _, exists := g.nodes[prefix]; !exists {
		if atomic.LoadUint64(&g.nodeCount) >= atomic.LoadUint64(&g.maxNodes) {
			return ErrNodeLimitExceeded
		}
		g.nodes[prefix] = &Node{
			Label:    prefix,
			Hash:     kmer.HashDouble(prefix),
			OutEdges: make(map[string]struct{}, EdgeBucketInitial),
			InEdges:  make(map[string]struct{}, EdgeBucketInitial),
		}
		atomic.AddUint64(&g.nodeCount, 1)
	}
	if _, exists := g.nodes[suffix]; !exists {
		if atomic.LoadUint64(&g.nodeCount) >= atomic.LoadUint64(&g.maxNodes) {
			return ErrNodeLimitExceeded
		}
		g.nodes[suffix] = &Node{
			Label:    suffix,
			Hash:     kmer.HashDouble(suffix),
			OutEdges: make(map[string]struct{}, EdgeBucketInitial),
			InEdges:  make(map[string]struct{}, EdgeBucketInitial),
		}
		atomic.AddUint64(&g.nodeCount, 1)
	}

	prefixNode := g.nodes[prefix]
	suffixNode := g.nodes[suffix]

	if _, exists := prefixNode.OutEdges[suffix]; !exists {
		if atomic.LoadUint64(&g.edgeCount) >= atomic.LoadUint64(&g.maxEdges) {
			return ErrEdgeLimitExceeded
		}
		prefixNode.OutEdges[suffix] = struct{}{}
		suffixNode.InEdges[prefix] = struct{}{}
		atomic.AddUint64(&g.edgeCount, 1)
	}

	atomic.AddUint32(&prefixNode.RefCount, 1)
	atomic.AddUint32(&suffixNode.RefCount, 1)

	return nil
}

func (g *Graph) AddSequence(seq string) error {
	kmers := kmer.Extract(seq, g.k)
	for _, k := range kmers {
		if err := g.AddKmer(k); err != nil {
			return err
		}
	}
	return nil
}

func (g *Graph) GetNode(key string) (*Node, bool) {
	g.mu.RLock()
	defer g.mu.RUnlock()
	n, ok := g.nodes[key]
	return n, ok
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

func (g *Graph) ForEachNode(fn func(label string, node *Node) bool) {
	g.mu.RLock()
	defer g.mu.RUnlock()
	for k, v := range g.nodes {
		if !fn(k, v) {
			return
		}
	}
}

func (g *Graph) Snapshot() *GraphSnapshot {
	g.mu.RLock()
	defer g.mu.RUnlock()
	snap := &GraphSnapshot{
		Nodes: make(map[string]*NodeSnapshot, len(g.nodes)),
		K:     g.k,
	}
	for k, v := range g.nodes {
		outs := make([]string, 0, len(v.OutEdges))
		for e := range v.OutEdges {
			outs = append(outs, e)
		}
		ins := make([]string, 0, len(v.InEdges))
		for e := range v.InEdges {
			ins = append(ins, e)
		}
		snap.Nodes[k] = &NodeSnapshot{
			Label:    v.Label,
			Hash:     v.Hash,
			OutEdges: outs,
			InEdges:  ins,
		}
	}
	return snap
}

type NodeSnapshot struct {
	Label    string
	Hash     kmer.DoubleHash
	OutEdges []string
	InEdges  []string
}

type GraphSnapshot struct {
	Nodes map[string]*NodeSnapshot
	K     int
}

func (g *Graph) VerifyIntegrity() int {
	g.mu.RLock()
	defer g.mu.RUnlock()
	collisions := 0
	hashToLabels := make(map[kmer.DoubleHash][]string)
	for label, node := range g.nodes {
		if !node.Hash.Equals(kmer.HashDouble(label)) {
			log.Printf("INTEGRITY ERROR: Node %q has mismatched hash (stored=%+v, recomputed=%+v)",
				label, node.Hash, kmer.HashDouble(label))
			collisions++
		}
		hashToLabels[node.Hash] = append(hashToLabels[node.Hash], label)
	}
	for h, labels := range hashToLabels {
		if len(labels) > 1 {
			log.Printf("COLLISION WARNING: DoubleHash %+v shared by %d labels: %v", h, len(labels), labels)
			collisions += len(labels) - 1
		}
	}
	return collisions
}
