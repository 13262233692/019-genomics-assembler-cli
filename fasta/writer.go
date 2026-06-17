package fasta

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

type Sequence struct {
	Header   string
	Sequence string
}

type Writer struct {
	writer io.Writer
	buf    *bufio.Writer
}

func NewWriter(w io.Writer) *Writer {
	return &Writer{writer: w, buf: bufio.NewWriter(w)}
}

func NewFileWriter(path string) (*Writer, error) {
	f, err := os.Create(path)
	if err != nil {
		return nil, err
	}
	return NewWriter(f), nil
}

func (w *Writer) Write(seq *Sequence) error {
	if _, err := fmt.Fprintf(w.buf, ">%s\n", seq.Header); err != nil {
		return err
	}
	lineWidth := 80
	seqLen := len(seq.Sequence)
	for i := 0; i < seqLen; i += lineWidth {
		end := i + lineWidth
		if end > seqLen {
			end = seqLen
		}
		if _, err := w.buf.WriteString(seq.Sequence[i:end] + "\n"); err != nil {
			return err
		}
	}
	return nil
}

func (w *Writer) Flush() error {
	return w.buf.Flush()
}
