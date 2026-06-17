package cmap

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

type EnzymeSite struct {
	Position    float64
	StdDev      float64
	Coverage    int
	Occurrence  int
	SiteID      int
	LabelChan   int
}

type ContigMap struct {
	ID           string
	Length       float64
	NumSites     int
	Sites        []EnzymeSite
	SitePositions []float64
}

type CMAP struct {
	FilePath    string
	EnzymeName  string
	Version     string
	ContigMaps  map[string]*ContigMap
	OrderedIDs  []string
}

func OpenCMAP(path string) (*CMAP, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ParseCMAP(f)
}

func ParseCMAP(r io.Reader) (*CMAP, error) {
	cmap := &CMAP{
		ContigMaps: make(map[string]*ContigMap),
		OrderedIDs: make([]string, 0),
	}

	scanner := bufio.NewScanner(r)
	scanner.Buffer(make([]byte, 1024*1024), 100*1024*1024)

	headerCols := make(map[string]int)
	dataStart := false

	for scanner.Scan() {
		line := scanner.Text()
		line = strings.TrimSpace(line)

		if line == "" {
			continue
		}

		if strings.HasPrefix(line, "#") {
			if strings.Contains(line, "enzyme name") || strings.Contains(line, "EnzymeName") {
				fields := strings.Fields(line)
				for i, f := range fields {
					if strings.Contains(f, "enzyme") || strings.Contains(f, "Enzyme") {
						if i+1 < len(fields) {
							cmap.EnzymeName = strings.Trim(fields[i+1], "\"':,")
						}
					}
				}
			}
			if strings.Contains(line, "CMAP version") || strings.Contains(line, "format version") {
				fields := strings.Fields(line)
				for _, f := range fields {
					if strings.HasPrefix(f, "v") || strings.Contains(f, ".") {
						cmap.Version = f
					}
				}
			}
			if strings.HasPrefix(line, "# ") && strings.Contains(line, "CMapId") {
				headerLine := strings.TrimPrefix(line, "# ")
				cols := strings.Split(headerLine, "\t")
				if len(cols) < 3 {
					cols = strings.Fields(headerLine)
				}
				for i, col := range cols {
					col = strings.TrimSpace(col)
					headerCols[col] = i
				}
				dataStart = true
			}
			continue
		}

		if !dataStart {
			continue
		}

		fields := strings.Split(line, "\t")
		if len(fields) < 3 {
			fields = strings.Fields(line)
		}
		if len(fields) < 3 {
			continue
		}

		cmapID := getFieldDef(fields, headerCols, 0, "CMapId", "CompntId")
		if cmapID == "" {
			continue
		}

		contigLen, _ := strconv.ParseFloat(getFieldDef(fields, headerCols, 1, "ContigLength", "Length"), 64)
		numSites, _ := strconv.Atoi(getFieldDef(fields, headerCols, 2, "NumSites", "Sites"))

		siteID, _ := strconv.Atoi(getFieldDef(fields, headerCols, 3, "SiteID", "Id"))
		labelChan, _ := strconv.Atoi(getFieldDef(fields, headerCols, 4, "LabelChannel", "Channel"))
		position, _ := strconv.ParseFloat(getFieldDef(fields, headerCols, 5, "Position", "Pos"), 64)
		stdDev, _ := strconv.ParseFloat(getFieldDef(fields, headerCols, 6, "StdDev", "SD"), 64)
		coverage, _ := strconv.Atoi(getFieldDef(fields, headerCols, 7, "Coverage", "Cov"))
		occurrence, _ := strconv.Atoi(getFieldDef(fields, headerCols, 8, "Occurrence", "Occ"))

		cm, exists := cmap.ContigMaps[cmapID]
		if !exists {
			cm = &ContigMap{
				ID:            cmapID,
				Length:        contigLen,
				NumSites:      numSites,
				Sites:         make([]EnzymeSite, 0, numSites),
				SitePositions: make([]float64, 0, numSites),
			}
			cmap.ContigMaps[cmapID] = cm
			cmap.OrderedIDs = append(cmap.OrderedIDs, cmapID)
		}

		site := EnzymeSite{
			Position:   position,
			StdDev:     stdDev,
			Coverage:   coverage,
			Occurrence: occurrence,
			SiteID:     siteID,
			LabelChan:  labelChan,
		}
		cm.Sites = append(cm.Sites, site)
		cm.SitePositions = append(cm.SitePositions, position)
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading CMAP: %v", err)
	}

	for _, cm := range cmap.ContigMaps {
		cm.NumSites = len(cm.Sites)
	}

	return cmap, nil
}

func getField(fields []string, headerCols map[string]int, names ...string) string {
	for _, name := range names {
		if idx, ok := headerCols[name]; ok && idx < len(fields) {
			return strings.TrimSpace(fields[idx])
		}
	}
	for _, name := range names {
		for h, idx := range headerCols {
			if strings.EqualFold(h, name) && idx < len(fields) {
				return strings.TrimSpace(fields[idx])
			}
		}
	}
	return ""
}

func getFieldDef(fields []string, headerCols map[string]int, defaultIdx int, names ...string) string {
	val := getField(fields, headerCols, names...)
	if val != "" {
		return val
	}
	if defaultIdx >= 0 && defaultIdx < len(fields) {
		return strings.TrimSpace(fields[defaultIdx])
	}
	return ""
}

func (c *CMAP) ContigCount() int {
	return len(c.ContigMaps)
}

func (c *CMAP) TotalSites() int {
	total := 0
	for _, cm := range c.ContigMaps {
		total += len(cm.Sites)
	}
	return total
}

func (c *CMAP) GetContig(id string) (*ContigMap, bool) {
	cm, ok := c.ContigMaps[id]
	return cm, ok
}

func (c *CMAP) LongestContig() *ContigMap {
	var longest *ContigMap
	for _, id := range c.OrderedIDs {
		cm := c.ContigMaps[id]
		if longest == nil || cm.Length > longest.Length {
			longest = cm
		}
	}
	return longest
}
