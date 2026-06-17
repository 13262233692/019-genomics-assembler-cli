package kmer

import (
	"hash/fnv"
	"strings"
)

const validBases = "ACGTacgt"

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
