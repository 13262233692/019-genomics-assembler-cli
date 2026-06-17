package debruijn

import (
	"testing"
)

func TestNewGraph(t *testing.T) {
	g := NewGraph(5)
	if g.K() != 5 {
		t.Errorf("Expected k=5, got %d", g.K())
	}
	if g.NodeCount() != 0 {
		t.Errorf("Expected 0 nodes, got %d", g.NodeCount())
	}
	if g.EdgeCount() != 0 {
		t.Errorf("Expected 0 edges, got %d", g.EdgeCount())
	}
}

func TestAddKmer(t *testing.T) {
	g := NewGraph(5)
	err := g.AddKmer("ABCDE")
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}

	if g.NodeCount() != 2 {
		t.Errorf("Expected 2 nodes, got %d", g.NodeCount())
	}
	if g.EdgeCount() != 1 {
		t.Errorf("Expected 1 edge, got %d", g.EdgeCount())
	}

	if g.OutDegree("ABCD") != 1 {
		t.Errorf("Expected out-degree 1 for ABCD, got %d", g.OutDegree("ABCD"))
	}
	if g.InDegree("BCDE") != 1 {
		t.Errorf("Expected in-degree 1 for BCDE, got %d", g.InDegree("BCDE"))
	}
}

func TestAddSequence(t *testing.T) {
	g := NewGraph(4)
	err := g.AddSequence("ACGTGACTAG")
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}

	expectedNodes := 8
	if g.NodeCount() != expectedNodes {
		t.Errorf("Expected %d nodes, got %d", expectedNodes, g.NodeCount())
	}

	startNodes := 0
	endNodes := 0
	g.ForEachNode(func(label string, node *Node) bool {
		if len(node.InEdges) == 0 {
			startNodes++
		}
		if len(node.OutEdges) == 0 {
			endNodes++
		}
		return true
	})
	if startNodes != 1 {
		t.Errorf("Expected 1 start node, got %d", startNodes)
	}
	if endNodes != 1 {
		t.Errorf("Expected 1 end node, got %d", endNodes)
	}
}

func TestDuplicateKmers(t *testing.T) {
	g := NewGraph(4)
	_ = g.AddKmer("ACGT")
	_ = g.AddKmer("ACGT")
	_ = g.AddKmer("ACGT")

	if g.NodeCount() != 2 {
		t.Errorf("Expected 2 nodes after duplicate kmers, got %d", g.NodeCount())
	}
	if g.EdgeCount() != 1 {
		t.Errorf("Expected 1 edge (not 3), got %d", g.EdgeCount())
	}
	if g.OutDegree("ACG") != 1 {
		t.Errorf("Expected out-degree 1 (not 3), got %d", g.OutDegree("ACG"))
	}
}

func TestGetNode(t *testing.T) {
	g := NewGraph(3)
	_ = g.AddKmer("ABC")

	node, ok := g.GetNode("AB")
	if !ok {
		t.Error("Expected to find node AB")
	}
	if len(node.OutEdges) != 1 {
		if _, exists := node.OutEdges["BC"]; !exists {
			t.Errorf("Expected out-edge to BC")
		}
	}

	_, ok = g.GetNode("XX")
	if ok {
		t.Error("Expected not to find node XX")
	}
}

func TestNodeLimits(t *testing.T) {
	g := NewGraphWithLimits(4, 2, 100)
	err := g.AddKmer("ACGT")
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	if g.NodeCount() < 2 {
		t.Errorf("Expected at least 2 nodes, got %d", g.NodeCount())
	}

	err = g.AddKmer("TTTT")
	if err != ErrNodeLimitExceeded {
		t.Errorf("Expected ErrNodeLimitExceeded, got %v", err)
	}
}

func TestEdgeLimits(t *testing.T) {
	g := NewGraphWithLimits(4, 100, 1)
	_ = g.AddKmer("AAAA")
	err := g.AddKmer("AAAT")
	if err != ErrEdgeLimitExceeded {
		t.Errorf("Expected ErrEdgeLimitExceeded, got %v", err)
	}
}

func TestVerifyIntegrity(t *testing.T) {
	g := NewGraph(4)
	sequences := []string{
		"ACGTACGTAGCT",
		"GCTAGCATCGAT",
		"TTTAAACCCGGG",
	}
	for _, s := range sequences {
		_ = g.AddSequence(s)
	}
	collisions := g.VerifyIntegrity()
	if collisions != 0 {
		t.Errorf("Expected 0 collisions, got %d", collisions)
	}
}

func TestHighlyRepetitiveSequence(t *testing.T) {
	g := NewGraph(5)
	rep := "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	err := g.AddSequence(rep)
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	if g.NodeCount() != 1 {
		t.Errorf("Expected 1 node for all-A repetitive sequence, got %d", g.NodeCount())
	}
	if g.EdgeCount() != 1 {
		t.Errorf("Expected 1 self-referencing edge, got %d", g.EdgeCount())
	}
	node, ok := g.GetNode("AAAA")
	if !ok {
		t.Fatal("Expected node AAAA")
	}
	if _, selfLoop := node.OutEdges["AAAA"]; !selfLoop {
		t.Error("Expected self-loop on AAAA node")
	}
}

func TestSnapshot(t *testing.T) {
	g := NewGraph(4)
	_ = g.AddSequence("ACGTGACTAG")
	snap := g.Snapshot()
	if snap.K != g.K() {
		t.Errorf("Snapshot K mismatch: %d vs %d", snap.K, g.K())
	}
	if len(snap.Nodes) != g.NodeCount() {
		t.Errorf("Snapshot node count mismatch: %d vs %d", len(snap.Nodes), g.NodeCount())
	}
	for label, ns := range snap.Nodes {
		node, ok := g.GetNode(label)
		if !ok {
			t.Errorf("Snapshot contains node %q not in graph", label)
			continue
		}
		if ns.Label != node.Label {
			t.Errorf("Node %q label mismatch in snapshot", label)
		}
		if len(ns.OutEdges) != len(node.OutEdges) {
			t.Errorf("Node %q out-edge count mismatch", label)
		}
	}
}
