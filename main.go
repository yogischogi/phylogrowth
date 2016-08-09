// Package phylogrowth calculates population growth from a phylogenetic tree.
package main

import (
	"bytes"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"time"

	"github.com/yogischogi/phylogrowth/histogram"
	"github.com/yogischogi/phylogrowth/phylotree"
)

func main() {
	// Command line flags.
	var (
		treein   = flag.String("treein", "", "Input filename for phylogenetic tree (.txt).")
		treeout  = flag.String("treeout", "", "Output filename for phylogenetic tree in TXT format.")
		csvout   = flag.String("csvout", "", "Output filename for histogram data in CSV format.")
		txtout   = flag.String("txtout", "", "Output filename for histogram data in TXT format.")
		pngout   = flag.String("pngout", "", "Output filename for PNG image. Needs Gnuplot to be installed.")
		step     = flag.Int("step", 100, "Step length of histogram intervals in years.")
		subclade = flag.String("subclade", "", "Selects a specific branch of the tree.")
	)
	flag.Parse()

	// Load phylogenetic tree from file.
	if *treein == "" {
		fmt.Printf("No filename for input tree specified.\r\n")
		os.Exit(1)
	}
	tree, err := phylotree.NewFromFile(*treein)
	if err != nil {
		fmt.Printf("Error reading tree from file, %v.\r\n", err)
		os.Exit(1)
	}

	// Select subclade.
	if *subclade != "" {
		tree = tree.Subclade(*subclade)
		if tree == nil {
			fmt.Printf("Error, could not find specified subclade %s.\r\n", *subclade)
			os.Exit(1)
		}
	}

	// Save resulting tree to file or print it out.
	if *treeout != "" {
		date := time.Now().Format("2006 Jan 2")
		var buffer bytes.Buffer
		buffer.WriteString("// This tree was created by the phylogrowth program: https://github.com/yogischogi/phylogrowth\r\n")
		buffer.WriteString("// Command used:\r\n// ")
		for _, arg := range os.Args {
			buffer.WriteString(arg)
			buffer.WriteString(" ")
		}
		buffer.WriteString("\r\n")
		buffer.WriteString("// " + date + "\r\n\r\n")
		buffer.WriteString(tree.String())
		err := ioutil.WriteFile(*treeout, buffer.Bytes(), os.ModePerm)
		if err != nil {
			fmt.Printf("Error writing tree to file, %v.\r\n", err)
		}
	}

	// Write histogram data to file.
	histo := histogram.New(tree.TMRCAs(), *step)
	if *csvout != "" {
		err := histo.WriteCSV(*csvout)
		if err != nil {
			fmt.Printf("Error writing histogram in CSV format to file, %v.\r\n", err)
		}
	}
	if *txtout != "" {
		err := histo.WriteTXT(*txtout)
		if err != nil {
			fmt.Printf("Error writing histogram in TXT format to file, %v.\r\n", err)
		}
	}
	if *pngout != "" {
		err := histo.WritePNG(*pngout, *subclade)
		if err != nil {
			fmt.Printf("Error writing histogram as PNG image format to file, %v.\r\n", err)
		}
	}
}
