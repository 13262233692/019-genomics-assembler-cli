package fastq

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strings"
)

type Read struct {
	Name     string
	Sequence string
	Quality  string
}

type Parser struct {
	scanner *bufio.Scanner
}

func NewParser(r io.Reader) *Parser {
	return &Parser{scanner: bufio.NewScanner(r)}
}

func OpenReader(path string) (io.ReadCloser, error) {
	if path == "" || path == "-" {
		return os.Stdin, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	if strings.HasSuffix(strings.ToLower(path), ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			f.Close()
			return nil, err
		}
		return &gzipReadCloser{gz, f}, nil
	}
	return f, nil
}

type gzipReadCloser struct {
	*gzip.Reader
	file *os.File
}

func (g *gzipReadCloser) Close() error {
	err1 := g.Reader.Close()
	err2 := g.file.Close()
	if err1 != nil {
		return err1
	}
	return err2
}

func (p *Parser) Next() (*Read, error) {
	name, err := p.readLine()
	if err != nil {
		if err == io.EOF {
			return nil, nil
		}
		return nil, err
	}
	if name == "" {
		return nil, nil
	}
	name = strings.TrimPrefix(name, "@")
	seq, err := p.readLine()
	if err != nil {
		return nil, err
	}
	plus, err := p.readLine()
	if err != nil {
		return nil, err
	}
	_ = plus
	qual, err := p.readLine()
	if err != nil {
		return nil, err
	}
	return &Read{Name: name, Sequence: seq, Quality: qual}, nil
}

func (p *Parser) readLine() (string, error) {
	if p.scanner.Scan() {
		return strings.TrimSpace(p.scanner.Text()), nil
	}
	if err := p.scanner.Err(); err != nil {
		return "", err
	}
	return "", io.EOF
}
