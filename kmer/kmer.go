package kmer

import (
	"hash/fnv"
	"strings"
)

const validBases = "ACGTacgt"

type DoubleHash struct {
	H1 uint64
	H2 uint64
}

func IsValid(seq string) bool {
	for _, b := range seq {
		if !strings.ContainsRune(validBases, b) {
			return false
		}
	}
	return true
}

func Extract(sequence string, k int) []string {
	if k <= 0 || len(sequence) < k {
		return nil
	}
	count := len(sequence) - k + 1
	kmers := make([]string, 0, count)
	for i := 0; i <= len(sequence)-k; i++ {
		kmer := sequence[i : i+k]
		if IsValid(kmer) {
			kmers = append(kmers, strings.ToUpper(kmer))
		}
	}
	return kmers
}

func Hash(kmer string) uint64 {
	h := fnv.New64a()
	h.Write([]byte(kmer))
	return h.Sum64()
}

func Hash2(kmer string) uint64 {
	data := []byte(kmer)
	var x uint64 = 14695981039346656037
	for _, b := range data {
		x ^= uint64(b)
		x *= 1099511628211
	}
	return x
}

func HashDouble(kmer string) DoubleHash {
	return DoubleHash{
		H1: Hash(kmer),
		H2: Hash2(kmer),
	}
}

func (d DoubleHash) Equals(other DoubleHash) bool {
	return d.H1 == other.H1 && d.H2 == other.H2
}

func Canonical(kmer string) string {
	rev := ReverseComplement(kmer)
	if kmer < rev {
		return kmer
	}
	return rev
}

func ReverseComplement(seq string) string {
	complement := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
		'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
	}
	result := make([]byte, len(seq))
	for i := 0; i < len(seq); i++ {
		result[len(seq)-1-i] = complement[seq[i]]
	}
	return string(result)
}

func EncodeKmer(kmer string) uint64 {
	var val uint64
	for i := 0; i < len(kmer) && i < 32; i++ {
		val <<= 2
		switch kmer[i] {
		case 'A', 'a':
			val |= 0
		case 'C', 'c':
			val |= 1
		case 'G', 'g':
			val |= 2
		case 'T', 't':
			val |= 3
		}
	}
	return val
}

func DecodeKmer(val uint64, k int) string {
	bases := []byte("ACGT")
	result := make([]byte, k)
	for i := k - 1; i >= 0; i-- {
		result[i] = bases[val&3]
		val >>= 2
	}
	return string(result)
}
