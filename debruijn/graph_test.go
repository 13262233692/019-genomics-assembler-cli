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
}

func TestAddKmer(t *testing.T) {
	g := NewGraph(5)
	g.AddKmer("ABCDE")

	if g.NodeCount() != 2 {
		t.Errorf("Expected 2 nodes, got %d", g.NodeCount())
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
	g.AddSequence("ACGTGACTAG")

	expectedNodes := 8
	if g.NodeCount() != expectedNodes {
		t.Errorf("Expected %d nodes, got %d", expectedNodes, g.NodeCount())
	}

	startNodes := 0
	endNodes := 0
	nodes := g.AllNodes()
	for _, node := range nodes {
		if len(node.InEdges) == 0 {
			startNodes++
		}
		if len(node.OutEdges) == 0 {
			endNodes++
		}
	}
	if startNodes != 1 {
		t.Errorf("Expected 1 start node, got %d", startNodes)
	}
	if endNodes != 1 {
		t.Errorf("Expected 1 end node, got %d", endNodes)
	}
}

func TestDuplicateKmers(t *testing.T) {
	g := NewGraph(4)
	g.AddKmer("ACGT")
	g.AddKmer("ACGT")
	g.AddKmer("ACGT")

	if g.NodeCount() != 2 {
		t.Errorf("Expected 2 nodes after duplicate kmers, got %d", g.NodeCount())
	}
	if g.OutDegree("ACG") != 1 {
		t.Errorf("Expected out-degree 1 (not 3), got %d", g.OutDegree("ACG"))
	}
}

func TestGetNode(t *testing.T) {
	g := NewGraph(3)
	g.AddKmer("ABC")

	node, ok := g.GetNode("AB")
	if !ok {
		t.Error("Expected to find node AB")
	}
	if len(node.OutEdges) != 1 || node.OutEdges[0] != "BC" {
		t.Errorf("Expected out-edge to BC, got %v", node.OutEdges)
	}

	_, ok = g.GetNode("XX")
	if ok {
		t.Error("Expected not to find node XX")
	}
}
