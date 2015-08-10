# Copyright 2015 Dorival Pedroso and Raul Durand. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

all: GenVtu ConvertGofemMat MatTable PlotLrm LocCmDriver ResidPlot Msh2vtu
.PHONY: GenVtu ConvertGofemMat MatTable PlotLrm LocCmDriver ResidPlot Msh2vtu

ConvertGofemMat: ConvertGofemMat.go
	go build -o /tmp/gofem/ConvertGofemMat ConvertGofemMat.go && mv /tmp/gofem/ConvertGofemMat $(GOPATH)/bin/

GenVtu: GenVtu.go
	go build -o /tmp/gofem/GenVtu GenVtu.go && mv /tmp/gofem/GenVtu $(GOPATH)/bin/

MatTable: MatTable.go
	go build -o /tmp/gofem/MatTable MatTable.go && mv /tmp/gofem/MatTable $(GOPATH)/bin/

PlotLrm: PlotLrm.go
	go build -o /tmp/gofem/PlotLrm PlotLrm.go && mv /tmp/gofem/PlotLrm $(GOPATH)/bin/

LocCmDriver: LocCmDriver.go
	go build -o /tmp/gofem/LocCmDriver LocCmDriver.go && mv /tmp/gofem/LocCmDriver $(GOPATH)/bin/

ResidPlot: ResidPlot.go
	go build -o /tmp/gofem/ResidPlot ResidPlot.go && mv /tmp/gofem/ResidPlot $(GOPATH)/bin/

Msh2vtu: Msh2vtu.go
	go build -o /tmp/gofem/Msh2vtu Msh2vtu.go && mv /tmp/gofem/Msh2vtu $(GOPATH)/bin/
