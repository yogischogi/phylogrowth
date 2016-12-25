// Package phylotree implements a phylogenetic tree with TMRCA estimates.
package phylotree

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"
	"unicode"
)

// Sample is line of text that represents a genetic sample from a single person.
type Sample struct {
	// Text is the text line used in the text representation of the tree.
	Text string
	ID   string
}

// newSample creates a new Sample from a textual representation.
func newSample(text string) (Sample, error) {
	result := Sample{Text: text}
	IDidx := strings.Index(text, "id:")
	IDidx += 3
	if IDidx > len([]rune(text)) {
		return result, errors.New(fmt.Sprintf("could not parse sample: %s", text))
	}
	result.ID = extractID(text[IDidx:])
	return result, nil
}

func (s *Sample) String() string {
	return fmt.Sprintf("id:%s", s.ID)
}

// Clade represents a phylogenetic clade.
type Clade struct {
	// Text is the text line used in the text representation of the tree.
	Text  string
	TMRCA float64
	// SNP is the name of the first SNP that characterizes this clade.
	SNP       string
	Subclades []Clade
	Samples   []Sample
}

// newClade creates a new Clade from a textual representation.
func newClade(text string) (Clade, error) {
	result := Clade{Text: text}
	result.SNP = extractID(text)

	// Search for TMRCA.
	tmrcaIdx := strings.Index(text, "TMRCA")
	if tmrcaIdx >= 0 && len(text) > tmrcaIdx+6 {
		tmrcaPart := text[tmrcaIdx+6:]
		if len(tmrcaPart) > 0 {
			tokens := strings.Fields(tmrcaPart)
			if len(tokens) > 0 {
				tmrca, err := strconv.ParseFloat(tokens[0], 64)
				if err != nil {
					msg := fmt.Sprintf("could not convert TMRCA value to float: %s", tokens[0])
					return result, errors.New(msg)
				}
				result.TMRCA = tmrca
			}
		}
	}
	return result, nil
}

// TMRCAs returns the TMRCAs for all new lineages.
func (c *Clade) TMRCAs() []float64 {
	result := make([]float64, 0)
	if c.TMRCA != 0 {
		branches := len(c.Subclades) + len(c.Samples)
		for i := 0; i < branches-1; i++ {
			result = append(result, c.TMRCA)
		}
		for _, subclade := range c.Subclades {
			subclade.appendTMRCAs(&result)
		}
	}
	return result
}

// appendTMRCAs appends the TMRCAs for lineages and samples
// from this clade to the result.
func (c *Clade) appendTMRCAs(result *[]float64) {
	if c.TMRCA != 0 {
		branches := len(c.Subclades) + len(c.Samples)
		for i := 0; i < branches-1; i++ {
			*result = append(*result, c.TMRCA)
		}
		for _, subclade := range c.Subclades {
			subclade.appendTMRCAs(result)
		}
	}
}

func (c *Clade) AddSample(sample Sample) {
	if c.Samples == nil {
		c.Samples = make([]Sample, 0)
	}
	c.Samples = append(c.Samples, sample)
}

func (c *Clade) AddSubclade(clade Clade) {
	if c.Subclades == nil {
		c.Subclades = make([]Clade, 0)
	}
	c.Subclades = append(c.Subclades, clade)
}

// NewFromFile parses a text file to create a tree.
// The return value Clade is the root node of the tree.
func NewFromFile(filename string) (*Clade, error) {
	lines := make([]lineInfo, 0)

	// Open file
	infile, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer infile.Close()

	// Read lines
	lineNo := 0
	scanner := bufio.NewScanner(infile)
	for scanner.Scan() {
		lineNo++
		text := stripComments(scanner.Text())
		// Use only valid lines and remove new samples.
		if text != "" && strings.Index(text, "new") == -1 {
			indent := countSpaces(text)
			text = strings.TrimSpace(text)
			lines = append(lines, lineInfo{lineNo: lineNo, indent: indent, text: text})
		}
	}
	if scanner.Err() != nil {
		return nil, scanner.Err()
	}
	lines = unitizeSingleLineages(lines)
	if len(lines) == 0 {
		return nil, errors.New("no valid lines in file, nothing to do")
	}

	// Build tree by parsing lines.
	root, err := newClade(lines[0].text)
	if err != nil {
		return nil, errors.New(fmt.Sprintf("invalid root element, %s", err))
	}
	err = parseTree(&root, lines[0].indent, lines[1:])
	if err != nil {
		return nil, errors.New(fmt.Sprintf("parsing tree, %s", err))
	}
	return &root, nil
}

func (c *Clade) String() string {
	var buffer bytes.Buffer
	c.prettyPrint(&buffer, 0)
	return buffer.String()
}

// Subclade returns the subclade that contains searchTerm.
func (c *Clade) Subclade(cladeName string) *Clade {
	var result *Clade
	if c.contains(cladeName) {
		result = c
	} else {
		for i, _ := range c.Subclades {
			if c.Subclades[i].contains(cladeName) {
				result = &c.Subclades[i]
				break
			} else {
				result = c.Subclades[i].Subclade(cladeName)
				if result != nil {
					break
				}
			}
		}
	}
	return result
}

// contains checks if the text representation of this clade
// contains cladeName.
// The method only returns true if the found cladeName is
// terminated by a non digit. This should make sure that the
// found cladeName is not part of a longer SNP name.
func (c *Clade) contains(cladeName string) bool {
	idx := strings.Index(c.Text, cladeName)
	if idx >= 0 {
		// Check if the search term is terminated or part of a longer SNP name.
		termIdx := idx + len(cladeName)
		if termIdx < len(c.Text) && unicode.IsDigit(rune(c.Text[termIdx])) {
			// Search term is part of longer SNP name.
			return false
		} else {
			return true
		}
	}
	return false
}

// prettyPrint prints a formatted version of the clade c
// into buffer. indent is the indentation for the root node.
func (c *Clade) prettyPrint(buffer *bytes.Buffer, indent int) {
	// Write this Element.
	for i := 0; i < indent; i++ {
		buffer.WriteString("\t")
	}
	buffer.WriteString(fmt.Sprintf("%s, TMRCA: %.0f\r\n", c.SNP, c.TMRCA))

	// Write Samples.
	for _, sample := range c.Samples {
		for i := 0; i < indent+1; i++ {
			buffer.WriteString("\t")
		}
		buffer.WriteString(sample.String())
		buffer.WriteString("\r\n")
	}
	// Write Subclades.
	for _, clade := range c.Subclades {
		clade.prettyPrint(buffer, indent+1)
	}
}

// lineInfo is a helper struct for parsing a tree in text format.
type lineInfo struct {
	lineNo int
	indent int
	text   string
}

// parseTree parses a tree in text format with white space indentations.
// The function works recursively and adds all new subclades and samples
// to the parent clade. indent is the indentation of the parent clade
// in the text file.
func parseTree(parent *Clade, indent int, lines []lineInfo) error {
	childIndent := -1
	for i, _ := range lines {
		switch {
		case lines[i].indent <= indent:
			// Return if indentation shows beginning of next block.
			return nil
		case childIndent == -1 && lines[i].indent > indent:
			// Determine the indentation of the child block.
			childIndent = lines[i].indent
			fallthrough
		case lines[i].indent == childIndent:
			// Parse child elements.
			if strings.Contains(lines[i].text, "id:") {
				sample, err := newSample(lines[i].text)
				if err != nil {
					msg := fmt.Sprintf("line: %d, %s", lines[i].lineNo, err)
					return errors.New(msg)
				}
				parent.AddSample(sample)
			} else {
				// Child is Clade element.
				clade, err := newClade(lines[i].text)
				if err != nil {
					msg := fmt.Sprintf("line: %d, %s", lines[i].lineNo, err)
					return errors.New(msg)
				}
				parseTree(&clade, lines[i].indent, lines[i+1:])
				parent.AddSubclade(clade)
			}
		}
	}
	return nil
}

// unitizeSingleLineages converts lines that are below a
// starred line, for example CTS4528*, into normal lines.
func unitizeSingleLineages(lines []lineInfo) []lineInfo {
	result := make([]lineInfo, 0, len(lines))
	indent := -1
	for _, line := range lines {
		if line.text[len(line.text)-1:] == "*" {
			indent = line.indent
			continue
		}
		switch {
		case indent == -1:
			// standard case
			result = append(result, line)
		case line.indent > indent:
			// Line belongs to a branch that ends with *.
			line.indent = indent
			result = append(result, line)
		case line.indent <= indent:
			// Line is outside the starred block.
			indent = -1
			result = append(result, line)
		}
	}
	return result
}

// stripComments removes comments from a line of text.
// A comment starts with // like in many programming languages.
// If a line contains only a comment or only whitespace characters
// an empty string is returned.
func stripComments(line string) string {
	result := line
	// Remove comment from line.
	idxComment := strings.Index(line, "//")
	if idxComment >= 0 {
		result = line[0:idxComment]
	}
	// Check if line contains only whitespaces.
	isEmpty := true
	for _, c := range result {
		if unicode.IsSpace(c) == false {
			isEmpty = false
			break
		}
	}
	if isEmpty {
		return ""
	} else {
		return result
	}
}

// extractID extracts the first ID from a line of text.
// An ID starts with one or more letters and ends with one or more digits.
func extractID(text string) string {
	var buffer bytes.Buffer
	startLetters := false
	startDigits := false
	for _, c := range text {
		isDigit := unicode.IsDigit(c)
		isLetter := unicode.IsLetter(c)
		if startLetters == false && !isLetter {
			continue
		} else if isDigit {
			startDigits = true
			buffer.WriteRune(c)
		} else if isLetter && startDigits == false {
			startLetters = true
			buffer.WriteRune(c)
		} else {
			break
		}
	}
	return buffer.String()
}

// countSpaces counts the white spaces at the beginning
// of a line.
func countSpaces(line string) int {
	result := 0
	for _, c := range line {
		if unicode.IsSpace(c) {
			result++
		} else {
			break
		}
	}
	return result
}
