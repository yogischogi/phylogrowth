// Package histogram provides methods to convert a series of numbers into histogram data.
package histogram

import (
	"bytes"
	"errors"
	"fmt"
	"io/ioutil"
	"math"
	"os"
	"os/exec"
)

// Histogram holds the values for the x and y scale as well
// as the original values.
type Histogram struct {
	values []float64
	x      []int
	y      []int
}

// New creates a new Histogram from a series of values
// and the step length of an interval.
func New(values []float64, step int) *Histogram {
	h := Histogram{values: values}
	roundedValues := make([]int, len(values))
	for i, v := range values {
		roundedValues[i] = int(math.Floor(float64(v)/float64(step)+0.5)) * step
	}

	min, max := minMax(roundedValues)
	nxValues := int((max-min)/step) + 1

	// Calculate x values.
	h.x = make([]int, nxValues)
	for i := 0; i < nxValues; i++ {
		h.x[i] = min + i*step
	}

	// Calculate y values.
	h.y = make([]int, nxValues)
	for _, v := range roundedValues {
		h.y[(v-min)/step]++
	}
	return &h
}

func (h *Histogram) WriteCSV(filename string) error {
	var buffer bytes.Buffer
	for i, _ := range h.x {
		line := fmt.Sprintf("%d,%d\r\n", h.x[i], h.y[i])
		buffer.WriteString(line)
	}
	err := ioutil.WriteFile(filename, buffer.Bytes(), os.ModePerm)
	return err
}

// WriteTXT writes the x and the y values of the histogram
// into a text file. Tabs are used to separate the columns.
func (h *Histogram) WriteTXT(filename string) error {
	var buffer bytes.Buffer
	for i, _ := range h.x {
		line := fmt.Sprintf("%d\t%d\r\n", h.x[i], h.y[i])
		buffer.WriteString(line)
	}
	err := ioutil.WriteFile(filename, buffer.Bytes(), os.ModePerm)
	return err
}

// WritePNG creates a PNG image from the data and writes it to a file.
// The subclade name is used for the title of the image.
//
// This method uses Gnuplot to create the image.
func (h *Histogram) WritePNG(filename string, subclade string) error {
	// Write temporary data file.
	dataFile, err := ioutil.TempFile("", "tempdata")
	if err != nil {
		return err
	}
	defer os.Remove(dataFile.Name())
	err = h.WriteTXT(dataFile.Name())
	if err != nil {
		return err
	}

	// Create temporary script for Gnuplot.
	scriptFile, err := ioutil.TempFile("", "tempscript")
	if err != nil {
		return err
	}
	defer os.Remove(scriptFile.Name())
	script := "reset\n" +
		"set terminal png\n" +
		"set xlabel \"Years before present\"\n" +
		"set ylabel \"Population growth\"\n" +
		"set title \"" + subclade + " Population Growth\"\n" +
		"set key below\n" +
		"set grid\n" +
		"set xrange[] reverse\n" +
		"plot \"" + dataFile.Name() + "\" using 1:2 title \"\" with linespoints"
	err = ioutil.WriteFile(scriptFile.Name(), []byte(script), os.ModePerm)
	if err != nil {
		return err
	}

	// Invoke Gnuplot.
	gnuplot := exec.Command("gnuplot", scriptFile.Name())
	var img bytes.Buffer
	gnuplot.Stdout = &img
	err = gnuplot.Run()
	if err != nil {
		return errors.New("could not run Gnuplot, maybe it is not installed")
	}
	// Write image.
	err = ioutil.WriteFile(filename, img.Bytes(), os.ModePerm)
	if err != nil {
		return err
	}
	return nil
}

func (h *Histogram) String() string {
	var buffer bytes.Buffer
	for i, _ := range h.x {
		line := fmt.Sprintf("%d, %d\r\n", h.x[i], h.y[i])
		buffer.WriteString(line)
	}
	return buffer.String()
}

// minMax returns the minimum and maximum values of values.
func minMax(values []int) (min, max int) {
	min = values[0]
	max = values[0]
	for _, v := range values {
		if v < min {
			min = v
		}
		if v > max {
			max = v
		}
	}
	return min, max
}
