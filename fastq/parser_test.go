package fastq

import (
	"strings"
	"testing"
)

func TestParser(t *testing.T) {
	input := "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nAAAA\n"
	reader := strings.NewReader(input)
	parser := NewParser(reader)

	read1, err := parser.Next()
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	if read1 == nil {
		t.Fatal("Expected read1, got nil")
	}
	if read1.Name != "read1" {
		t.Errorf("Expected name 'read1', got '%s'", read1.Name)
	}
	if read1.Sequence != "ACGT" {
		t.Errorf("Expected sequence 'ACGT', got '%s'", read1.Sequence)
	}
	if read1.Quality != "IIII" {
		t.Errorf("Expected quality 'IIII', got '%s'", read1.Quality)
	}

	read2, err := parser.Next()
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	if read2 == nil {
		t.Fatal("Expected read2, got nil")
	}
	if read2.Name != "read2" {
		t.Errorf("Expected name 'read2', got '%s'", read2.Name)
	}

	read3, err := parser.Next()
	if err != nil {
		t.Fatalf("Unexpected error: %v", err)
	}
	if read3 != nil {
		t.Error("Expected nil read at end of file")
	}
}
