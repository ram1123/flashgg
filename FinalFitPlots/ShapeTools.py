from sys import stdout, stderr
import ROOT

from HiggsAnalysis.CombinedLimit.ModelTools import ModelBuilder

class ShapeBuilder(ModelBuilder):
    def __init__(self,datacard,options):
        ModelBuilder.__init__(self,datacard,options) 
        if not datacard.hasShapes: 
            raise RuntimeError, "You're using a ShapeBuilder for a model that has no shapes"
        if options.libs:
            for lib in options.libs:
                ROOT.gSystem.Load(lib)
    ## ------------------------------------------
    ## -------- ModelBuilder interface ----------
    ## ------------------------------------------
    	self.wsp = None
    def doObservables(self):
        if (self.options.verbose > 1): stderr.write("Using shapes: qui si parra' la tua nobilitate\n")
        self.prepareAllShapes();
        if len(self.DC.bins) > 1:
            strexpr="CMS_channel[" + ",".join(["%s=%d" % (l,i) for i,l in enumerate(self.DC.bins)]) + "]";
            self.doVar(strexpr);
            self.out.binCat = self.out.cat("CMS_channel");
            if self.options.verbose: stderr.write("Will use category 'CMS_channel' to identify the %d channels\n" % self.out.binCat.numTypes())
            self.out.obs = ROOT.RooArgSet()
            self.out.obs.add(self.out.binVars)
            self.out.obs.add(self.out.binCat)
        else:
            self.out.obs = self.out.binVars
        self.doSet("observables",self.out.obs)
        if len(self.DC.obs) != 0: 
            self.doCombinedDataset()
    def doIndividualModels(self):
        for b in self.DC.bins:
            pdfs   = ROOT.RooArgList(); bgpdfs   = ROOT.RooArgList()
            coeffs = ROOT.RooArgList(); bgcoeffs = ROOT.RooArgList()
            for p in self.DC.exp[b].keys(): # so that we get only self.DC.processes contributing to this bin
                (pdf,coeff) = (self.getPdf(b,p), self.out.function("n_exp_bin%s_proc_%s" % (b,p)))
                extranorm = self.getExtraNorm(b,p)
                if extranorm:
                    self.doObj("n_exp_final_bin%s_proc_%s" % (b,p), "prod", "n_exp_bin%s_proc_%s, %s" % (b,p, extranorm))
                    coeff = self.out.function("n_exp_final_bin%s_proc_%s" % (b,p))                    
                pdfs.add(pdf); coeffs.add(coeff)
                if not self.DC.isSignal[p]:
                    bgpdfs.add(pdf); bgcoeffs.add(coeff)
            sum_s = ROOT.RooAddPdf("pdf_bin%s"       % b, "",   pdfs,   coeffs)
            sum_b = ROOT.RooAddPdf("pdf_bin%s_bonly" % b, "", bgpdfs, bgcoeffs)
            if len(self.DC.systs):
                # rename the pdfs
                sum_s.SetName("pdf_bin%s_nuis" % b); sum_b.SetName("pdf_bin%s_bonly_nuis" % b)
                # now we multiply by all the nuisances, but avoiding nested products
                # so we first make a list of all nuisances plus the RooAddPdf
                sumPlusNuis_s = ROOT.RooArgList(self.out.nuisPdfs); sumPlusNuis_s.add(sum_s)
                sumPlusNuis_b = ROOT.RooArgList(self.out.nuisPdfs); sumPlusNuis_b.add(sum_b)
                # then make RooProdPdf and import it
                pdf_s = ROOT.RooProdPdf("pdf_bin%s"       % b, "", sumPlusNuis_s) 
                pdf_b = ROOT.RooProdPdf("pdf_bin%s_bonly" % b, "", sumPlusNuis_b) 
                self.out._import(pdf_s, ROOT.RooFit.RenameConflictNodes(b))
                self.out._import(pdf_b, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
            else:
                self.out._import(sum_s, ROOT.RooFit.RenameConflictNodes(b))
                self.out._import(sum_b, ROOT.RooFit.RecycleConflictNodes(), ROOT.RooFit.Silence())
    def doCombination(self):
        ## Contrary to Number-counting models, here each channel PDF already contains the nuisances
        ## So we just have to build the combined pdf
        if len(self.DC.bins) > 1:
            for (postfixIn,postfixOut) in [ ("","_s"), ("_bonly","_b") ]:
                simPdf = ROOT.RooSimultaneous("model"+postfixOut, "model"+postfixOut, self.out.binCat)
                for b in self.DC.bins:
                    simPdf.addPdf(self.out.pdf("pdf_bin%s%s" % (b,postfixIn)), b)
                self.out._import(simPdf)
        else:
            self.out._import(self.out.pdf("pdf_bin%s"       % self.DC.bins[0]).clone("model_s"), ROOT.RooFit.Silence())
            self.out._import(self.out.pdf("pdf_bin%s_bonly" % self.DC.bins[0]).clone("model_b"), ROOT.RooFit.Silence())
  
    ## --------------------------------------
    ## -------- High level helpers ----------
    ## --------------------------------------
    def prepareAllShapes(self):
        shapeTypes = []; shapeBins = []; shapeObs = {}
        for ib,b in enumerate(self.DC.bins):
            for p in [self.options.dataname]+self.DC.exp[b].keys():
                if len(self.DC.obs) == 0 and p == self.options.dataname: continue
                if p != self.options.dataname and self.DC.exp[b][p] == 0: continue
                shape = self.getShape(b,p); norm = 0;
                if shape == None: # counting experiment
                    if not self.out.var("CMS_fakeObs"): 
                        self.doVar("CMS_fakeObs[0,1]");
                        self.doSet("CMS_fakeObsSet","CMS_fakeObs");
                        shapeObs["CMS_fakeObsSet"] = self.out.set("CMS_fakeObsSet")
                    if p == self.options.dataname:
                        shapeTypes.append("RooDataHist")
                    else:
                        shapeTypes.append("RooAbsPdf");
                elif shape.ClassName().startswith("TH1"):
                    shapeTypes.append("TH1"); shapeBins.append(shape.GetNbinsX())
                    norm = shape.Integral()
                elif shape.InheritsFrom("RooDataHist"):
                    shapeTypes.append("RooDataHist"); 
                    shapeBins.append(shape.numEntries())
                    shapeObs[self.argSetToString(shape.get(0))] = shape.get(0)
                    norm = shape.sumEntries()
                elif shape.InheritsFrom("RooDataSet"):
                    shapeTypes.append("RooDataSet"); 
                    shapeObs[self.argSetToString(shape.get(0))] = shape.get(0)
                    norm = shape.sumEntries()
                elif shape.InheritsFrom("TTree"):
                    shapeTypes.append("TTree"); 
                elif shape.InheritsFrom("RooAbsPdf"):
                    shapeTypes.append("RooAbsPdf");
                else: raise RuntimeError, "Currently supporting only TH1s, RooDataHist and RooAbsPdfs"
                if norm != 0:
                    if p == self.options.dataname:
                        if len(self.DC.obs):
                            if self.DC.obs[b] == -1: self.DC.obs[b] = norm
                            elif abs(norm-self.DC.obs[b]) > 0.01:
                                raise RuntimeError, "Mismatch in normalizations for observed data in bin %s: text %f, shape %f" % (b,self.DC.obs[b],norm)
                    else:
                        if self.DC.exp[b][p] == -1: self.DC.exp[b][p] = norm
                        elif abs(norm-self.DC.exp[b][p]) > 0.01*max(1,self.DC.exp[b][p]): 
                            raise RuntimeError, "Mismatch in normalizations for bin %s, process %s: rate %f, shape %f" % (b,p,self.DC.exp[b][p],norm)
        if shapeTypes.count("TH1") == len(shapeTypes):
            self.out.allTH1s = True
            self.out.mode    = "binned"
            self.out.maxbins = max(shapeBins)
            if self.options.verbose: stderr.write("Will use binning variable 'x' with %d bins\n" % self.out.maxbins)
            self.doVar("x[0,%d]" % self.out.maxbins); self.out.var("x").setBins(self.out.maxbins)
            self.out.binVar = self.out.var("x")
            self.out.binVars = ROOT.RooArgSet(self.out.binVar)
        elif shapeTypes.count("RooDataSet") > 0 or shapeTypes.count("TTree") > 0:
            self.out.mode = "unbinned"
            if self.options.verbose: stderr.write("Will try to work with unbinned datasets\n")
            if self.options.verbose: stderr.write("Observables: %s\n" % str(shapeObs.keys()))
            if len(shapeObs.keys()) != 1:
                self.out.binVars = ROOT.RooArgSet()
                for obs in shapeObs.values():
                     self.out.binVars.add(obs, False)
            else:
                self.out.binVars = shapeObs.values()[0]
            self.out._import(self.out.binVars)
        else:
            self.out.mode = "binned"
            if self.options.verbose: stderr.write("Will try to make a binned dataset\n")
            if self.options.verbose: stderr.write("Observables: %s\n" % str(shapeObs.keys()))
            if len(shapeObs.keys()) != 1:
                raise RuntimeError, "There's more than once choice of observables: %s\n" % str(shapeObs.keys())
            self.out.binVars = shapeObs.values()[0]
            self.out._import(self.out.binVars)
    def doCombinedDataset(self):
        if len(self.DC.bins) == 1:
            data = self.getData(self.DC.bins[0],self.options.dataname).Clone(self.options.dataname)
            self.out._import(data)
            return
        if self.out.mode == "binned":
            combiner = ROOT.CombDataSetFactory(self.out.obs, self.out.binCat)
            for b in self.DC.bins: combiner.addSet(b, self.getData(b,self.options.dataname))
            self.out.data_obs = combiner.done(self.options.dataname,self.options.dataname)
            self.out._import(self.out.data_obs)
        elif self.out.mode == "unbinned":
            combiner = ROOT.CombDataSetFactory(self.out.obs, self.out.binCat)
            for b in self.DC.bins: combiner.addSet(b, self.getData(b,self.options.dataname))
            self.out.data_obs = combiner.doneUnbinned(self.options.dataname,self.options.dataname)
            self.out._import(self.out.data_obs)
        else: raise RuntimeException, "Only combined datasets are supported"
    ## -------------------------------------
    ## -------- Low level helpers ----------
    ## -------------------------------------
    def getShape(self,channel,process,syst="",_fileCache={},_cache={}):
        if _cache.has_key((channel,process,syst)): 
            if self.options.verbose: print "recyling (%s,%s,%s) -> %s\n" % (channel,process,syst,_cache[(channel,process,syst)].GetName())
            return _cache[(channel,process,syst)];
        bentry = None
        if self.DC.shapeMap.has_key(channel): bentry = self.DC.shapeMap[channel]
        elif self.DC.shapeMap.has_key("*"):   bentry = self.DC.shapeMap["*"]
        else: raise KeyError, "Shape map has no entry for channel '%s'" % (channel)
        names = []
        if bentry.has_key(process): names = bentry[process]
        elif bentry.has_key("*"):   names = bentry["*"]
        elif self.DC.shapeMap["*"].has_key(process): names = self.DC.shapeMap["*"][process]
        elif self.DC.shapeMap["*"].has_key("*"):     names = self.DC.shapeMap["*"]["*"]
        else: raise KeyError, "Shape map has no entry for process '%s', channel '%s'" % (process,channel)
        if len(names) == 1 and names[0] == "FAKE": return None
        if syst != "": names = [names[0], names[2]]
        else:          names = [names[0], names[1]]
        strmass = "%d" % self.options.mass if self.options.mass % 1 == 0 else str(self.options.mass)
        finalNames = [ x.replace("$PROCESS",process).replace("$CHANNEL",channel).replace("$SYSTEMATIC",syst).replace("$MASS",strmass) for x in names ]
        if not _fileCache.has_key(finalNames[0]): _fileCache[finalNames[0]] = ROOT.TFile.Open(finalNames[0])
        file = _fileCache[finalNames[0]]; objname = finalNames[1]
        if not file: raise RuntimeError, "Cannot open file %s (from pattern %s)" % (finalNames[0],names[0])
        if ":" in objname: # workspace:obj or ttree:xvar or th1::xvar
	    print "Going for workspace"
            (wname, oname) = objname.split(":")
	    if self.wsp is None: self.wsp = file.Get(wname)
	    #self.wsp = file.Get(wname)
            if not self.wsp: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (wname,finalNames[0],names[1],names[0])
            if self.wsp.ClassName() == "RooWorkspace":
                ret = self.wsp.data(oname)
                if not ret: ret = self.wsp.pdf(oname)
                if not ret: raise RuntimeError, "Object %s in workspace %s in file %s does not exist or it's neither a data nor a pdf" % (oname, wname, finalNames[0])
                ret.SetName("shape_%s_%s%s" % (process,channel, "_"+syst if syst else ""))
                _cache[(channel,process,syst)] = ret
                if not syst:
                  normname = "%s_norm" % (oname)
                  norm =self. wsp.arg(normname)
                  if norm: 
                    norm.SetName("shape_%s_%s%s_norm" % (process,channel, "_"))
                    self.out._import(norm, ROOT.RooFit.RecycleConflictNodes()) 
                if self.options.verbose: print "import (%s,%s) -> %s\n" % (finalNames[0],objname,ret.GetName())
                return ret;
            elif self.wsp.ClassName() == "TTree":
                ##If it is a tree we will convert it in RooDataSet . Then we can decide if we want to build a
                ##RooKeysPdf or if we want to use it as an unbinned dataset 
                if not self.wsp: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (wname,finalNames[0],names[1],names[0])
                self.doVar("%s[%f,%f]" % (oname,self.wsp.GetMinimum(oname),self.wsp.GetMaximum(oname)))
                #Check if it is weighted
                self.doVar("__WEIGHT__[0.,1000.]")
                rds = ROOT.RooDataSet("shape_%s_%s%s" % (process,channel, "_"+syst if syst else ""), "shape_%s_%s%s" % (process,channel, "_"+syst if syst else ""),self.wsp,ROOT.RooArgSet(self.out.var(oname)),"","__WEIGHT__")
                rds.var = oname
                _cache[(channel,process,syst)] = rds
                if self.options.verbose: print "import (%s,%s) -> %s\n" % (finalNames[0],wname,rds.GetName())
                return rds
            elif self.wsp.InheritsFrom("TH1"):
                ##If it is a Histogram we will convert it in RooDataSet preserving the bins 
                if not self.wsp: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (wname,finalNames[0],names[1],names[0])
                name = "shape_%s_%s%s" % (process,channel, "_"+syst if syst else "")
                # don't make it twice
                for X in _neverDelete:
                    if X.InheritsFrom("TNamed") and X.GetName() == name: return X
                self.doVar("%s[%f,%f]" % (oname,self.wsp.GetXaxis().GetXmin(),self.wsp.GetXaxis().GetXmax()))
                rds = ROOT.RooDataHist(name, name, ROOT.RooArgList(self.out.var(oname)), self.wsp)
                rds.var = oname
                if self.options.verbose: stderr.write("import (%s,%s) -> %s\n" % (finalNames[0],wname,rds.GetName()))
                _neverDelete.append(rds)
                return rds
            else:
                raise RuntimeError, "Object %s in file %s has unrecognized type %s" (wname, finalNames[0], self.wsp.ClassName())
        else: # histogram
            ret = file.Get(objname);
            if not ret: raise RuntimeError, "Failed to find %s in file %s (from pattern %s, %s)" % (objname,finalNames[0],names[1],names[0])
            ret.SetName("shape_%s_%s%s" % (process,channel, "_"+syst if syst else ""))
            if self.options.verbose: print "import (%s,%s) -> %s\n" % (finalNames[0],objname,ret.GetName())
            _cache[(channel,process,syst)] = ret
            return ret
    def getData(self,channel,process,syst="",_cache={}):
        return self.shape2Data(self.getShape(channel,process,syst),channel,process)
    def getPdf(self,channel,process,_cache={}):
        if _cache.has_key((channel,process)): return _cache[(channel,process)]
        shapeNominal = self.getShape(channel,process)
        nominalPdf = self.shape2Pdf(shapeNominal,channel,process)
        if shapeNominal == None: return nominalPdf # no point morphing a fake shape
        morphs = []; shapeAlgo = None
        for (syst,pdf,args,errline) in self.DC.systs:
            if not "shape" in pdf: continue
            if shapeAlgo != None and pdf != shapeAlgo: raise RuntimeError, "You can use only one morphing algorithm for a given shape"
            shapeAlgo = pdf
            if errline[channel][process] != 0:
                shapeUp   = self.getShape(channel,process,syst+"Up")
                shapeDown = self.getShape(channel,process,syst+"Down")
                if shapeUp.ClassName()   != shapeNominal.ClassName(): raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst" % (channel,process,syst)
                if shapeDown.ClassName() != shapeNominal.ClassName(): raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst" % (channel,process,syst)
                morphs.append((syst,errline[channel][process],self.shape2Pdf(shapeUp,channel,process),self.shape2Pdf(shapeDown,channel,process)))
        if len(morphs) == 0: return nominalPdf
        if shapeAlgo == "shapeN": stderr.write("Warning: the shapeN implementation in RooStats and L&S are different\n")
        pdfs = ROOT.RooArgList(nominalPdf)
        coeffs = ROOT.RooArgList()
        minscale = 1
        for (syst,scale,pdfUp,pdfDown) in morphs:
            pdfs.add(pdfUp); pdfs.add(pdfDown);
            if scale == 1:
                coeffs.add(self.out.var(syst))
            else: # must scale it :-/
                coeffs.add(self.doObj("%s_scaled_%s_%s" % (syst,channel,process), "prod","%s, %s" % (scale,syst)))
                if scale < minscale: minscale = scale
        qrange = minscale; qalgo = 0;
        if shapeAlgo[-1] == "*": 
            qalgo = 100
            shapeAlgo = shapeAlgo[:-1]
        if shapeAlgo == "shapeL": qrange = 0;
        elif shapeAlgo == "shapeN": qalgo = -1;
        _cache[(channel,process)] = ROOT.VerticalInterpPdf("shape_%s_%s_morph" % (channel,process), "", pdfs, coeffs, qrange, qalgo)
        return _cache[(channel,process)]
    def getExtraNorm(self,channel,process):
        terms = []
        shapeNominal = self.getShape(channel,process)
        if shapeNominal == None: 
            # FIXME no extra norm for dummy pdfs (could be changed)
            return None
        if shapeNominal.InheritsFrom("RooAbsPdf"): 
            # return nominal multiplicative normalization constant
            normname = "shape_%s_%s%s_norm" % (process,channel, "_")
            if self.out.arg(normname): return normname
            else: return None
        normNominal = 0
        if shapeNominal.InheritsFrom("TH1"): normNominal = shapeNominal.Integral()
        elif shapeNominal.InheritsFrom("RooDataHist"): normNominal = shapeNominal.sumEntries()
        else: return None    
        for (syst,pdf,args,errline) in self.DC.systs:
            if "shape" not in pdf: continue
            if errline[channel][process] != 0:
                shapeUp   = self.getShape(channel,process,syst+"Up")
                shapeDown = self.getShape(channel,process,syst+"Down")
                if shapeUp.ClassName()   != shapeNominal.ClassName(): raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst" % (channel,process,syst)
                if shapeDown.ClassName() != shapeNominal.ClassName(): raise RuntimeError, "Mismatched shape types for channel %s, process %s, syst" % (channel,process,syst)
                kappaUp,kappaDown = 1,1
                if shapeNominal.InheritsFrom("TH1"):
                    kappaUp,kappaDown = shapeUp.Integral(),shapeDown.Integral()
                elif shapeNominal.InheritsFrom("RooDataHist"):
                    kappaUp,kappaDown = shapeUp.sumEntries(),shapeDown.sumEntries()
                kappaUp /=normNominal; kappaDown /= normNominal
                if abs(kappaUp-1) < 1e-3 and abs(kappaDown-1) < 1e-3: continue
                # if errline[channel][process] == <x> it means the gaussian should be scaled by <x> before doing pow
                # for convenience, we scale the kappas
                kappasScaled = [ pow(x, errline[channel][process]) for x in kappaDown,kappaUp ]
                terms.append("AsymPow(%f,%f,%s)" % (kappasScaled[0], kappasScaled[1], syst))
        return ",".join(terms) if terms else None;
    def shape2Data(self,shape,channel,process,_cache={}):
        if shape == None:
            name = "shape_%s_%s" % (channel,process)
            if not _cache.has_key(name):
                obs = ROOT.RooArgSet(self.out.var("CMS_fakeObs"))
                obs.setRealValue("CMS_fakeObs",0.5);
                if self.out.mode == "binned":
                    self.out.var("CMS_fakeObs").setBins(1)
                    rdh = ROOT.RooDataHist(name, name, obs)
                    rdh.set(obs, self.DC.obs[channel])
                    _cache[name] = rdh
                else:
                    rds = ROOT.RooDataSet(name, name, obs)
                    if self.DC.obs[channel] == float(int(self.DC.obs[channel])):
                        for i in range(int(self.DC.obs[channel])): rds.add(obs)
                    else:
                        rds.add(obs, self.DC.obs[channel])
                    _cache[name] = rds
            return _cache[name]
        if not _cache.has_key(shape.GetName()):
            if shape.ClassName().startswith("TH1"):
                rebinh1 = ROOT.TH1F(shape.GetName()+"_rebin", "", self.out.maxbins, 0.0, float(self.out.maxbins))
                for i in range(1,min(shape.GetNbinsX(),self.out.maxbins)+1): 
                    rebinh1.SetBinContent(i, shape.GetBinContent(i))
                rdh = ROOT.RooDataHist(shape.GetName(), shape.GetName(), ROOT.RooArgList(self.out.var("x")), rebinh1)
                self.out._import(rdh)
                _cache[shape.GetName()] = rdh
            elif shape.ClassName() in ["RooDataHist", "RooDataSet"]:
                return shape
            else: raise RuntimeError, "shape2Data not implemented for %s" % shape.ClassName()
        return _cache[shape.GetName()]
    def shape2Pdf(self,shape,channel,process,_cache={}):
        if shape == None:
            name = "shape_%s_%s" % (channel,process)
            if not _cache.has_key(name):
                _cache[name] = ROOT.RooUniform(name, name, ROOT.RooArgSet(self.out.var("CMS_fakeObs")))
            return _cache[name]
        if not _cache.has_key(shape.GetName()+"Pdf"):
            if shape.ClassName().startswith("TH1"):
                rdh = self.shape2Data(shape,channel,process)
                rhp = self.doObj("%sPdf" % shape.GetName(), "HistPdf", "{x}, %s" % shape.GetName())
                _cache[shape.GetName()+"Pdf"] = rhp
            elif shape.InheritsFrom("RooAbsPdf"):
                _cache[shape.GetName()+"Pdf"] = shape
            elif shape.InheritsFrom("RooDataHist"):
                rhp = ROOT.RooHistPdf("%sPdf" % shape.GetName(), "", self.out.binVars, shape) 
                self.out._import(rhp)
                _cache[shape.GetName()+"Pdf"] = rhp
            elif shape.InheritsFrom("RooDataSet"):
                rkp = ROOT.RooKeysPdf("%sPdf" % shape.GetName(), "", self.out.var(shape.var), shape,3,1.5); 
                self.out._import(rkp)
                _cache[shape.GetName()+"Pdf"] = rkp
            else: 
                raise RuntimeError, "shape2Pdf not implemented for %s" % shape.ClassName()
        return _cache[shape.GetName()+"Pdf"]
    def argSetToString(self,argset):
        names = []
        it = argset.createIterator()
        while True:
            arg = it.Next()
            if not arg: break
            names.append(arg.GetName())
        return ",".join(names)

