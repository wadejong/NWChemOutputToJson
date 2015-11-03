import json, sys, string, itertools
from collections import OrderedDict

class nwchemToJson:

  def __init__(self,processOrbitals='doOrbitals'):
    self.simulationEnv = {}
    self.simulationTime = {}
    self.calculations = []
    self.basis = basisObj()
    self.molecule = moleculeObj()

    self.calcSetup = {}
    self.calcRes = {}
    self.calcTask = {}

    self.setupCount = 0
    self.taskNumber = 0
    self.subTask = False

    self.doOrbitals = processOrbitals
   
  def convert(self,streamIn):
    mytasks = OrderedDict([
      ('echo of input',self.readInput),
      ('Job information',self.readEnv),
      ('Basis "ao',self.basis.readBasis),
      ('ECP       "e',self.basis.readEcp),
      ('Geometry "',self.molecule.readGeom),
      ('SCF Module',self.readScfDft),
      ('NWChem DFT Module',self.readScfDft),
      ('TDDFT Module',self.readTddft),
      ('Geometry Optimization',self.readGeomOpt),
      ('Frequency Analysis',self.readFreq),
      ('Property Module',self.readProp),
      ('Many-Electron Theory Module',self.readTce)
    ])

    collectingInput = False
    line = streamIn.readline()

    while line:
      if line.find('Input Module'): 
        collectingInput = True
      for myKey in mytasks.keys():
        if line.find(myKey)>=0: 
           myIndex = list(mytasks).index(myKey)
           if myIndex>4:
             self.calcTask = {}
             self.calcSetup = {}
             self.calcRes = {}
             self.taskNumber += 1
             self.setupCount += 1
             collectingInput = False
           if not collectingInput and myIndex in range (2,4): 
             break
           mytasks[myKey](line,streamIn)
      if line.find('Total times')>=0:
        self.simulationTime = {'simulationTime' : self.readTiming(line)}
        break
      line = streamIn.readline()

    return json.dumps({'simulation' : { 'simulationEnvironment' : self.simulationEnv,
                                        'calculations'          : self.calculations,
                                        'simulationTime'        : self.simulationTime }}, 
                                         indent = 2, separators=(',', ': '), ensure_ascii=False)

  def setMoleculeID(self):
    self.calcTask['id'] = 'calculation.'+str(self.taskNumber)
    self.calcTask['molecularFormula'] = self.molecule.molecularFormula
 
  def setSetup(self):
    if self.molecule.geomUpdated:
      self.calcSetup['molecule'] = self.molecule.molecule
      self.molecule.geomUpdated = False
    else:
      self.calcSetup['molecule'] = 'Molecule.'+str(self.molecule.molCount)
    if self.basis.basUpdated:
      self.calcSetup['basisSet'] = self.basis.basis
      self.basis.basUpdated = False
    else:
      self.calcSetup['basisSet'] = 'BasisSet.'+str(self.basis.basCount)
    self.calcSetup['id'] = 'calculationSetup.'+str(self.setupCount)

  def genFuncForMol(self):
    def generateFunctions(funcName,spherical):
      lValueList = {'s' : 0, 'p' : 1, 'd' : 2, 'f' : 3, 'g' : 4, 'h' : 5}
      lValue = lValueList[funcName[-1:]]
      functionList = [] 
      xyz = ['x','y','z']
      if spherical: 
        for i in range(-lValue, lValue+1): 
          functionList.append(str(i))
        return functionList
      def funcGen(lValue,funcName):
        if lValue == 0: 
          functionList.append(funcName)
        else:
          for comp in xyz:
            if xyz.index(comp) >= (xyz.index(funcName[-1:]) if len(funcName)>1 else -1): 
              funcGen(lValue-1,funcName+comp)
      funcGen(lValue,funcName)
      return functionList
    functionListForMolecule = []
    atomNum = 0
    for atom in self.molecule.molecule['atoms']:
      atomNum+=1
      for basis in self.basis.basis['basisFunctions']:
        if basis['elementLabel'] == atom['elementLabel'] and basis['basisSetType'] == 'orbitalBasis':
          spherical = basis['basisSetHarmonicType'] == 'spherical'
          for contr in basis['basisSetContraction']:
            atomString = str(atomNum)+' '+basis['elementLabel']+' '
            for x in generateFunctions(contr['basisSetShellType'],spherical): 
              functionListForMolecule.append(atomString+x)
    return functionListForMolecule

  def readInput(self,line,streamIn):
    inputData = ''
    line = streamIn.readline()
    while line:
      if line.find("====")>=0: 
        break
      inputData+=line.lstrip(' ').rstrip('\n')+';'
      line = streamIn.readline()
    self.simulationEnv['inputData'] = inputData.rstrip(';')

  def readEnv(self,line,streamIn):
    self.simulationEnv['programRun'] = "NWChem"
    envList = {
      'hostname'        : 'hostMachine',
      'date'            : 'runDate',
      'nwchem revision' : 'programVersion',
      'input'           : 'inputFileName',
      'nproc'           : 'processorCount'
    }
    while line:
      if line.find('time left')>=0: 
        break
      for envKey in envList.keys():
        if line.find(envKey)>=0: 
          self.simulationEnv[envList.get(envKey)] = line.split('=')[1].lstrip(' ').rstrip('\n')
      line = streamIn.readline()
  
  def readTaskTimes(self,line):
      vars = line.split()
      return { 'cpuTime' : float(vars[3].rstrip('s')), 'wallTime' : float(vars[5].rstrip('s')), 'units' : 'second' }
  
  def readTiming(self,line): 
    vars = line.split()
    return { 'cpuTime' : float(vars[3].rstrip('s')), 'wallTime' : float(vars[5].rstrip('s')), 'units' : 'second' }
        
  def readScfDft(self,line,streamIn):
    self.calcSetup['numberOfElectrons'] = 0
    self.calcSetup['molecularSpinMultiplicity'] = 1
    def closedShell(line):  
      self.calcSetup['numberOfElectrons'] = int(line.split()[3])*2
    def openShell(line): 
      val = line.split()[3]
      self.calcSetup['numberOfElectrons'] = int(val)+int(self.calcSetup['numberOfElectrons'])
      self.calcSetup['molecularSpinMultiplicity'] = int(val)/2+1
    def doCharge(line): 
      self.calcSetup['charge'] = int(float(line.split()[2]))
    def waveFuncSCF(line): 
      self.calcSetup['waveFunctionType'] = line.split()[2]
      self.calcSetup['waveFunctionTheory'] = 'Hartree-Fock'
    def waveFuncDFT(line): 
      if line.split()[2] == 'closed shell':
        self.calcSetup['waveFunctionType'] = 'RHF'
      else:
        self.calcSetup['waveFunctionType'] = 'UHF'
      self.calcSetup['waveFunctionTheory'] = 'Density Functional Theory'
    def inputVec(line): 
      self.calcSetup['inputVectors'] = line.split()[3]
    def outputVec(line): 
      self.calcSetup['outputVectors'] = line.split()[3]
    def alphaElec(line): 
      self.calcSetup['numberOfElectrons'] = int(line.split()[3]) + int(self.calcSetup['numberOfElectrons'])
    def betaElec(line): 
      val = line.split()[3]
      self.calcSetup['molecularSpinMultiplicity'] = abs(int(val)-int(self.calcSetup['numberOfElectrons']))/2+1
      self.calcSetup['numberOfElectrons'] = int(val) + int(self.calcSetup['numberOfElectrons'])
    def totalEn(line): 
      self.calcRes['totalEnergy'] = { 'value' : float(line.split()[4].replace('D','E')), 'units' : 'Hartree'}
    def oneEn(line): 
      self.calcRes['oneElectronEnergy'] = { 'value' : float(line.split()[3].replace('D','E')), 'units' : 'Hartree'}
    def oneEn2(line): 
      self.calcRes['oneElectronEnergy'] = { 'value' : float(line.split()[4].replace('D','E')), 'units' : 'Hartree'}
    def twoEn(line): 
      self.calcRes['twoElectronEnergy'] = { 'value' : float(line.split()[3].replace('D','E')), 'units' : 'Hartree'}
    def nucEn(line): 
      self.calcRes['nuclearRepulsionEnergy'] = { 'value' : float(line.split()[4].replace('D','E')), 'units' : 'Hartree'}
    def coulEn(line): 
      self.calcRes['coulombEnergy'] = {'value' : float(line.split()[3].replace('D','E')), 'units' : 'Hartree'}
    def xcEn(line): 
      self.calcRes['exchangeCorrelationEnergy'] = {'value' : float(line.split()[3].replace('D','E')), 'units' : 'Hartree'}
    def xEn(line): 
      self.calcRes['exchangeEnergy'] = {'value' : float(line.split()[3].replace('D','E')), 'units' : 'Hartree'}
    def cEn(line): 
      self.calcRes['correlationEnergy'] = {'value' : float(line.split()[3].replace('D','E')), 'units' : 'Hartree'}
    def s2Val(line): 
      self.calcRes['s2ExpectationValue'] = {'value' : float(line.split()[2]), 'units' : 'none'}
    def szVal(line): 
      self.calcRes['szExpectationValue'] = {'value' : float(line.split()[2]), 'units' : 'none'}
    def xcFunc(streamIn):
      exchange = []
      for _ in range(2): 
        line = streamIn.readline()
      while line.find('Grid Info')<0 and len(line.split())>1:
        if line.find('Method')>=0:
          exchange.append({'xcName' : line.lstrip(' ').rstrip('\n')})
        elif line.find('Correlation')>=0:
          vars = line.lstrip(' ').rstrip('\n').partition('Functional')
          fac = vars[2].split()[0]
          loc = ''
          if len(vars[2].split())>1 : 
            loc = vars[2].split()[1]
            exchange.append({'correlationTermName' : vars[0]+vars[1],
                             'correlationTermFactor' : float(fac),
                             'correlationTermLocality' : loc
                             })
          else:
            exchange.append({'correlationTermName' : vars[0]+vars[1],
                             'correlationTermFactor' : float(fac)
                             })
        else:
          if line.find('Exact')>=0: vars = line.lstrip(' ').rstrip('\n').partition('Exchange')
          else: vars = line.lstrip(' ').rstrip('\n').partition('Functional')
          fac = vars[2].split()[0]
          loc = ''
          if len(vars[2].split())>1 : 
            loc = vars[2].split()[1]
            exchange.append({'exchangeTermName' : vars[0]+vars[1],
                             'exchangeTermFactor' : float(fac),
                             'exchangeTermLocality' : loc
                            })
          else:
            exchange.append({'exchangeTermName' : vars[0]+vars[1],
                             'exchangeTermFactor' : float(fac)
                            })
        line = streamIn.readline()
      self.calcSetup['exchangeCorrelationFunctional'] = exchange
    def readOrbitals(line,streamIn):
      if self.doOrbitals == 'noOrbitals':
        return
      if line.find('Alpha')>=0:
        alphaBeta = True
        loopRange = 2
      else:
        alphaBeta = False
        loopRange = 1
      molecularOrbitals = {}
      molecularOrbitals['id'] = 'orbitalsMolecule.'+str(self.molecule.molCount)
      molecularOrbitals['atomicOrbitalDescriptions'] = self.genFuncForMol()
      molecularOrbital = []
      for i in range(loopRange):
        if loopRange>1:
          if i == 0: 
            spinLabel = '-alpha'
          else:
            spinLabel = '-beta'
        else:
          spinLabel = ''
        while line:
          if line.find('Vector')>=0:
            break
          line = streamIn.readline()
        while line:
          if line.find('Vector')>=0:
            orbital = {}
            coefficients = [0.0] * len(molecularOrbitals['atomicOrbitalDescriptions'])
            vars = line.replace('=',' ').split()
            orbital['id'] = 'molecularOrbital'+spinLabel+'.'+str(vars[1])+'.Mol.'+str(self.molecule.molCount)
            orbital['orbitalEnergy'] = { 'value' : float(vars[5].replace('D','E')), 'units' : 'Hartree' }
            orbital['orbitalOccupancy'] = float(vars[3].replace('D','E'))
            if line.find('Symmetry')>=0:
              orbital['orbitalSymmetry'] = vars[7]
            for _ in range(3):
              line = streamIn.readline()
            while True:
              vars = streamIn.readline().split()
              if len(vars) < 1:
                 break
              coefficients[int(vars[0])-1] = float(vars[1])
              if len(vars) > 6:
                coefficients[int(vars[5])-1] = float(vars[6])
            orbital['moCoefficients'] = coefficients
            molecularOrbital.append(orbital)
          else:
            break
          line = streamIn.readline()
      molecularOrbitals['molecularOrbital'] = molecularOrbital
      self.calcRes['molecularOrbitals'] = molecularOrbitals
    if not self.subTask:
      self.setMoleculeID()
      self.setSetup()
      self.calcTask['calculationType'] = 'energyCalculation'
    scfInp = {
      'closed shells'     : closedShell,
      'open shells'       : openShell,
      'Charge           :': doCharge,
      'charge          =' : doCharge,
      'wavefunction'      : waveFuncSCF,
      'Wavefunction type' : waveFuncDFT,
      'input vectors'     : inputVec,
      'output vectors'    : outputVec,
      'alpha elec'        : alphaElec,
      'beta elec'         : betaElec,
      'Alpha electrons'   : alphaElec,
      'Beta electrons'    : betaElec,
      'Total DFT ener'    : totalEn,
      'Total SCF ener'    : totalEn,
      'One-electron e'    : oneEn,
      'One electron e'    : oneEn2,
      'Two-electron e'    : twoEn,
      'Coulomb energy'    : coulEn,
      'Exchange-Corr'     : xcEn,
      'Exchange ener'     : xEn,
      'Correlation ener'  : cEn,
      'Nuclear repuls'    : nucEn,
      'S^2 ='             : s2Val,
      '<S2> ='            : s2Val,
      'Sz ='              : szVal
    }
    line = streamIn.readline()
    while line:
      if  line.find('Task  times')>=0 and not self.subTask:
        self.calcTask['calculationTime'] = self.readTaskTimes(line)
        break
      elif line.find('is already converged')>=0:
        line = streamIn.readline()
        self.calcRes['totalEnergy'] = { 'value' : streamIn.readline().split('=')[1], 'units' : 'Hartree'}
        break
      elif (line.find('Module')>=0 or line.find('Parallel integral file')>=0 or line.find('Line search')>=0 or line.find('Saving state')>=0) and self.subTask:
        break
      elif line.find('XC Information')>=0: 
        xcFunc(streamIn)
      elif line.find('Molecular Orbital Analysis')>=0: 
        readOrbitals(line,streamIn)
      else:
        for scfKey in scfInp.keys(): 
          if line.find(scfKey)>=0: 
            scfInp[scfKey](line)
      line = streamIn.readline()
    self.calcTask['calculationResults'] = self.calcRes
    if not self.subTask:
      self.calcTask['calculationSetup'] = self.calcSetup
      self.calculations.append(self.calcTask)

  def readTddft(self,line,streamIn):
    tddftStates = []
    self.setMoleculeID()
    self.setSetup()
    self.calcTask['calculationType'] = 'dftExcitedStates'
    while line.find('TDDFT Info')<0: line = streamIn.readline()
    line = streamIn.readline()
    self.calcTask['calculationType'] = streamIn.readline().split(':')[1].rstrip('\n')
    self.calcTask['waveFunctionType'] = streamIn.readline().split(':')[1].lstrip(' ').rstrip('\n')
    unrestricted = self.calcTask['waveFunctionType'].startswith('Unrestricted')
    while line.find('Ground state')<0: line = streamIn.readline()
    vars = line.split()
    self.calcRes['groundState'] = {'groundStateEnergy' : { 'value' : float(vars[3]), 'units' : 'Hartree'}, 'groundStateSymmetry' : vars[2]}
    if unrestricted: 
      self.calcRes['groundState'].update({'s2ExpectationValue': {'value' : float(line.split()[2]), 'units' : 'none'}})
    singlets = False
    triplets = False
    if not unrestricted : 
      singlets = self.calcTask['waveFunctionType'].find('singlets')>=0
      triplets = self.calcTask['waveFunctionType'].find('triplets')>=0
    if singlets or unrestricted:
      exstate = ''
      if singlets: exstate = '-singlet'
      while line:
        if line.find('Root  ')>=0: 
          vars = line.split()
          state = {'excitedState' : vars[1]+exstate, 'excitedStateSymmetry' : vars[3], 'excitationEnergy' : { 'value' : float(vars[4]), 'units' : 'Hartree'}}
          if unrestricted: state.update({'s2ExpectationValue':{'value':float(streamIn.readline().split()[2]), 'units':'none'}})
          for _ in range(4): 
            line = streamIn.readline()
          state.update({'dipoleOscillatorStrength' : { 'value' : float(streamIn.readline().split()[3]), 'units' : 'none'}})
          tddftStates.append(state)
        elif line.find('Target root')>=0:
          line = streamIn.readline()
          break
        line = streamIn.readline()
    if triplets: 
      exstate = '-triplet'
      while line:
        if line.find('Root')>=0: 
          vars = line.split()
          state = {'excitedState' : vars[1]+exstate, 'excitedStateSymmetry' : vars[3], 'excitationEnergy' : { 'value' : float(vars[4]), 'units' : 'Hartree'}}
          state.update({'dipoleOscillatorStrength' : { 'value' : 0.0, 'units' : 'none'}})
          tddftStates.append(state)
        elif line.find('Target root')>=0:
          line = streamIn.readline()
          break
        line = streamIn.readline()
    while line:
      if  line.find('Task  times')>=0: self.calcTask['calculationTime'] = self.readTaskTimes(line)
      line = streamIn.readline()
    self.calcRes['excitedStates'] = tddftStates
    self.calcTask['calculationResults'] = self.calcRes
    self.calculations.append(self.calcTask)
  
  def readGeomOpt(self,line,streamIn):
    self.subTask = True
    self.setMoleculeID()
    self.setSetup()
    self.calcTask['calculationSetup'] = self.calcSetup
    for _ in range(29):
        line = streamIn.readline()
    if line.find('Transition')>=0:
      self.calcTask['calculationType'] = 'geometrySaddlePoint'
    else:
      self.calcTask['calculationType'] = 'geometryOptimization'
    while line:
      if  line.find('Optimization converged')>=0:
        line = streamIn.readline() 
        while line:
          if line.find('Geometry "')>=0:
            self.molecule.readGeom(line,streamIn)
            self.molecule.geomUpdated = False
            break
          line=streamIn.readline()
        self.calcRes['molecule'] = self.molecule.molecule
      elif  line.find('Task  times')>=0:
        self.calcTask['calculationTime'] = self.readTaskTimes(line)
        break
      elif line.find('Module')>=0 and not line.find('Gradient')>=0: 
        self.readScfDft(line,streamIn)
      line = streamIn.readline()
    self.calcTask['calculationResults'] = self.calcRes
    self.calculations.append(self.calcTask)
    self.subTask = False
  
  def readFreq(self,line,streamIn):
    self.subTask = True
    self.setMoleculeID()
    self.setSetup()
    self.calcTask['calculationType'] = 'vibrationalModes'
    freq = []
    mods = [] 
    inten = []
    def tempFreq(line,streamIn): 
      self.calcRes['temperatureVibrations'] = { 'value' : float(line.split()[2].rstrip('K')), 'units' : 'Kelvin'}
    def freqScal(line,streamIn): 
      self.calcRes['frequencyScalingFactor'] = float(line.split()[4])
    def zpeBlock(line,streamIn):
      self.calcRes['zeroPointEnergyCorrection'] = { 'value' : float(line.split()[8]), 'units' : 'Hartree'}
      self.calcRes['thermalEnergyCorrection'] = { 'value' : float(streamIn.readline().split()[8]), 'units' : 'Hartree'}
      self.calcRes['thermalEnthalpyCorrection'] = { 'value' : float(streamIn.readline().split()[8]), 'units' : 'Hartree'}
    def entBlock(line,streamIn):
      entropy = {}
      vars = line.split()
      entropy['totalEntropy'] = { 'value' : float(vars[3]), 'units' : vars[4] }
      vars = streamIn.readline().split()
      entropy['translationalEntropyContribution'] = { 'value' : float(vars[3]), 'units' : vars[4], 'molecularWeight' :  vars[8].rstrip(')') }
      vars = streamIn.readline().split()
      entropy['rotationalEntropyContribution'] = { 'value' : float(vars[3]), 'units' : vars[4], 'symmetryNumber' :  vars[8].rstrip(')') }
      vars = streamIn.readline().split()
      entropy['vibrationalEntropyContribution'] = { 'value' : float(vars[3]), 'units' : vars[4] }
      self.calcRes['entropy'] = entropy
    def eigenVec(line,streamIn):
      for _ in range(2): 
        line = streamIn.readline()
      while line:
        numVal = len(line.split())
        if numVal == 0: 
          break
        line = streamIn.readline()
        vars = streamIn.readline().split()
        del vars[0]
        for x in vars:
          freq.append(x)
        for _ in range(2): 
          line = streamIn.readline()
        while line:
          if line.find('--')>=0 or len(line.split())<=1: 
            break
          else:
            vars = line.split()
            del vars[0]
            for x in vars:  
              mods.append(x)
          line = streamIn.readline()
        line = streamIn.readline()
    def intenVal(line,streamIn):
      for _ in range(2): 
        line = streamIn.readline()
      for _ in range(0,len(freq)):
        inten.append(streamIn.readline().split()[4])
    freqInp = {
      'Temperature'          : tempFreq,
      'frequency scaling'    : freqScal,
      'Zero-Point correct'   : zpeBlock,
      'Total Entropy'        : entBlock,
      'Projected Freq'       : eigenVec,
      'Projected Infra'      : intenVal
    }
    line = streamIn.readline()
    energyRead = True
    while line:
      if  line.find('Task  times')>=0:
        self.calcTask['calculationTime'] = self.readTaskTimes(line)
        break
      elif line.find('Module')>=0 and not line.find('CPHF')>=0: 
        if energyRead:
          self.readScfDft(line,streamIn)
          energyRead = False
      for freqKey in freqInp.keys():
        if line.find(freqKey)>=0:
          freqInp[freqKey](line,streamIn)
      line = streamIn.readline()
    normalModes = []
    numberOfModes = len(freq)
    minMod, maxMod = 0, 0
    while maxMod < (numberOfModes-1):
      minMod = maxMod
      maxMod = min(minMod+6,numberOfModes-1)
      if maxMod == 0: 
        break
      rangeMod = maxMod - minMod + 1
      for i in range (minMod, maxMod):
        modVec = []
        offset = minMod * numberOfModes
        if abs(float(freq[i])) > 0.0:
          for j in range(offset, offset + numberOfModes * rangeMod, rangeMod):
            modVec.append(float(mods[j]))
          normalModes.append({ 'id' : 'normalMode.'+str(i+1) , 
                               'normalModeFrequency' : { 'value' : float(freq[i]), 'units' : 'cm-1'}, 
                               'normalModeInfraRedIntensity' : { 'value' : float(inten[i]), 'units' : '(debye/angs)**2'},
                               'normalModeVector' : { 'value' :  modVec, 'units' : 'none'}
                             })
    self.calcRes['normalModes'] = normalModes
    self.calcTask['calculationSetup'] = self.calcSetup
    self.calcTask['calculationResults'] = self.calcRes
    self.calculations.append(self.calcTask)
    self.subTask = False
  
  def readProp(self,line,streamIn):
    self.subTask = True 
    self.setMoleculeID()
    self.setSetup()
    self.calcTask['calculationType'] = 'molecularProperties'
    def propDip(streamIn):
      myList = ['momentX', 'momentY', 'momentZ']
      dipProp = {}
      for _ in range(3): 
        line = streamIn.readline()
      vars = streamIn.readline().split()
      dipProp['expansionPoint'] = { 'value' : [float(vars[2]),float(vars[5]),float(vars[8])], 'units' : 'atomic units' }
      line = streamIn.readline()
      dipProp['totalMoment'] = { 'value' : float(streamIn.readline().split()[2]), 'units' : 'atomic units' }
      for i in range(3): 
        dipProp[myList[i]] = { 'value' : float(streamIn.readline().split()[1]), 'units' : 'atomic units' }
      molProps[0].update({'dipoleMoment' : dipProp})
      return line
    def propQuad(streamIn):
      myList = ['momentXX', 'momentYY', 'momentZZ', 'momentXY', 'momentXZ', 'momentYZ']
      quadProp = {}
      for _ in range(3): 
        line = streamIn.readline()
      vars = streamIn.readline().split()
      quadProp['expansionPoint'] = { 'value' : [float(vars[2]),float(vars[5]),float(vars[8])], 'units' : 'atomic units' }
      line = streamIn.readline()
      quadProp['diamagneticSusceptibility'] = { 'value' : float(streamIn.readline().split()[4]), 'units' : 'atomic units' }
      for _ in range(6): 
        line = streamIn.readline()
      for i in range(6): 
        quadProp[myList[i]] = { 'value' : float(streamIn.readline().split()[3]), 'units' : 'atomic units' }
      molProps[0].update({'quadrupoleMoment' : quadProp})
      return line
    def propESP(streamIn):
      for _ in range(4): 
        line = streamIn.readline()
      for atomNum in range(1,self.molecule.atomCount+1):
        vars = streamIn.readline().split()
        molProps[atomNum].update({'electrostaticPotential' : { 'value' : float(vars[5]), 'units' : 'atomic units' }})
        molProps[atomNum].update({'diamagneticShielding' : { 'value' : float(vars[6]), 'units' : 'atomic units' }})
      return line
    def propEfield(streamIn):
      myList = [ 'electricFieldX', 'electricFieldY', 'electricFieldZ', 'electricField' ]
      for _ in range(7): 
        line = streamIn.readline()
      for atomNum in range(1,self.molecule.atomCount+1):
        vars = streamIn.readline().split()
        for j in range(0,4):
          molProps[atomNum].update({myList[j] : { 'value' : float(vars[j+5]), 'units' : 'atomic units' }})
      return line
    def propEFG(streamIn):
      myIndex = [0,3,4,3,1,5,4,5,2]
      for _ in range(3): 
        line = streamIn.readline()
      for atomNum in range(1,self.molecule.atomCount+1):
        for _ in range(11): 
          line = streamIn.readline()
        vars = streamIn.readline().split()
        value = []
        efgProp = {}
        for j in range(9): 
          value.append(float(vars[myIndex[j]]))
        efgProp['efgTensor'] = { 'tensorValues' : value, 'units' : 'atomic units'}
        for _ in range(4): 
          line = streamIn.readline()
        vars = streamIn.readline().split()
        efgProp['efgPrincipalComponents'] = { 'values' : [float(vars[0]),float(vars[1]),float(vars[2])], 'units' : 'atomic units'}
        efgProp['efgAsymmetry'] = float(vars[3])
        line = streamIn.readline()
        value = []
        for i in range (3):
          vars = streamIn.readline().split()
          for j in range(3):
            value.insert(i*(j+1)+j,float(vars[j]))
        efgProp['efgProjectionVectors'] =  { 'vectorValues' : value, 'units' : 'none' }
        molProps[atomNum].update({'electricFieldGradient' : efgProp})
      return line
    def propShield(streamIn):
      line = streamIn.readline()
      chemshield = {}
      csNames = [ 'diamagneticComponentShieldingTensor' , 'paramagneticComponentShieldingTensor', 'shieldingTensor' ]
      linesSkipped = 0
      while line.find('Atom:')<0: 
        line = streamIn.readline()
      while line:
        if line.find('Atom:')>=0:
          atomNum = line.split()[1]
          line = streamIn.readline()
          for a in range(3):
            value = []
            for _ in range(3): 
              for x in streamIn.readline().split(): value.append(float(x))
            chemshield[csNames[a]] = { 'tensorValues' : value, 'units' : 'ppm' }
            for _ in range(2): line = streamIn.readline()
          chemshield['isotropicShielding'] = { 'value' : float(line.split()[2]), 'units' : 'ppm' }
          chemshield['shieldingAnisotropy'] = { 'value' : float(streamIn.readline().split()[2]), 'units' : 'ppm' }
          for _ in range(3): 
            line = streamIn.readline()
          value = streamIn.readline().split()
          chemshield['shieldingPrincipalComponents'] = { 'values' : [float(value[0]),float(value[1]),float(value[2])], 'units' : 'ppm' }
          line = streamIn.readline()
          value = []
          for i in range (0,3):
            vars = streamIn.readline().split()
            for j in range(0,3):
              value.insert(i*(j+1)+j,float(vars[j+1]))
          chemshield['shieldingProjectionVectors'] =  { 'vectorValues' : value, 'units' : 'none' }
          molProps[int(atomNum)].update({'electricFieldGradient' : chemshield})
          linesSkipped = 0
        elif linesSkipped>3: 
          break
        else:  
          linesSkipped += 1
        line = streamIn.readline()
      return line
    def propSpinSpin(streamIn):
      myList = [ 'Fermi Contact Term', 'Spin-Dipole Term', 'Fermi Contact - Spin-Dipole Cross Term',
                 'Paramagnetic Spin-Orbit Term', 'Diamagnetic Spin-Orbit Term', 'Spin-Spin Coupling Tensor']
      spinspin = []
      atomList = []
      line = streamIn.readline()
      while line.find('Indirect Spin-Spin Tensors')<0: 
        line=streamIn.readline()
      linesSkipped = 0
      while line:
        if line.find('Atom')>=0:
          vars = line.split()
          atom1 = vars[1].rstrip(':')
          atom2 = vars[5].rstrip(':')
          if atom1 not in atomList:
            atomList.append(atom1)
            spinspin.append({ 'atomicWeight' : float(vars[2].split('-')[0]), 'spinSpinCouplingPairs' : [] })
          if atom2 not in atomList:
            atomList.append(atom2)
            spinspin.append({ 'atomicWeight' : float(vars[6].split('-')[0]), 'spinSpinCouplingPairs' : [] })
          coupledAtom1 = { 'atom' : 'Atom.'+str(atom2)+'.Mol.'+str(molCount), 'atomicWeight' : float(vars[6].split('-')[0]) }
          coupledAtom2 = { 'atom' : 'Atom.'+str(atom1)+'.Mol.'+str(molCount), 'atomicWeight' : float(vars[2].split('-')[0]) }
          line=streamIn.readline()
          vars = streamIn.readline().split()
          if atom1 not in atomList: 
            spinspin.append({ 'nuclearGFactor' : float(vars[3])})
          if atom2 not in atomList: 
            spinspin.append({ 'nuclearGFactor' : float(vars[5])})
          coupledAtom1.update( {'nuclearGFactor' : float(vars[5])} )
          coupledAtom2.update( {'nuclearGFactor' : float(vars[3])} )
          for _ in range(3): 
            line = streamIn.readline()
          spinpair1 = {}
          for i in range(6):
            value = []
            for _ in range(3): 
              for x in streamIn.readline().split(): 
                value.append(float(x))
            if i<5 : 
              iso = streamIn.readline().split()[2]
            spinpair1[myList[i]] = {'tensorValues' : value, 'isotropicValue' : float(iso), 'units' : 'Hertz' } 
            for _ in range(2): 
              line = streamIn.readline()
          spinpair1['Isotropic Spin-Spin Coupling'] = {'value' : float(streamIn.readline().split()[4]), 'units' : 'Hertz'}
          spinpair2 = spinpair1
          spinpair1['coupledAtom'] = coupledAtom1
          spinpair2['coupledAtom'] = coupledAtom2
          spinspin[atomList.index(atom1)].update( {'spinSpinCouplingPairs' : spinpair1} )
          spinspin[atomList.index(atom2)].update( {'spinSpinCouplingPairs' : spinpair2} )
          linesSkipped = 0
        elif linesSkipped>3 and line.find('---')<0: 
          break
        else: 
          linesSkipped += 1
        for i in range(len(atomList)): 
          molProps[int(atomList[i])].update({'spinSpinCoupling' : spinspin[i-1]})
        line=streamIn.readline()
      return line
    propInp = {
      'Dipole Moment'       : propDip,
      'Quadrupole Moment'   : propQuad,
      'Electrostatic pot'   : propESP,
      'Electric field   '   : propEfield,
      'Electric field grad' : propEFG,
      'Chemical Shielding'  : propShield,
      'Indirect Spin-Spin'  : propSpinSpin
    }
    molProps = []
    molProps.append({ 'Molecule' : 'Molecule.'+str(self.molecule.molCount) })
    for i in range(1, self.molecule.atomCount+1):
      atomProp = {}
      atomProp['atom'] = 'Atom.'+str(i)+'.Mol.'+str(self.molecule.molCount)
      molProps.append(atomProp)
    line = streamIn.readline()
    while line:
      if line.find('Module')>=0 and not line.find('CPHF')>=0: 
         self.readScfDft(line,streamIn)
      for propKey in propInp.keys():
        if line.find(propKey)>=0: 
           line = propInp[propKey](streamIn)
      if  line.find('Task  times')>=0:
        self.calcTask['calculationTime'] = self.readTaskTimes(line)
        break
      line = streamIn.readline()
    for i in range(self.molecule.atomCount,-1,-1): 
      if len(molProps[i])<2: 
        del molProps[i]
    self.calcRes['molecularProperties'] = molProps
    self.calcTask['calculationSetup'] = self.calcSetup
    self.calcTask['calculationResults'] = self.calcRes
    self.calculations.append(self.calcTask)
    self.subTask = False

  def readTce(self,line,streamIn):
    self.setMoleculeID()
    self.setSetup()
    self.calcTask['calculationType'] = 'energyCalculation'
    while line.find('Number of processors')<0:
      line = streamIn.readline()
    for _ in range(2):
      line = streamIn.readline()
    if line.find('Restricted'):
      if line.find('open-shell'):
        self.calcSetup['waveFunctionType'] = 'ROHF'
      else:
        self.calcSetup['waveFunctionType'] = 'RHF'
    else:
      self.calcSetup['waveFunctionType'] = 'UHF'
    self.calcSetup['numberOfElectrons'] = int(streamIn.readline().split(':')[1])
    self.calcSetup['numberOfAlphaElectrons'] = int(streamIn.readline().split(':')[1])
    self.calcSetup['numberOfBetaElectrons'] = int(streamIn.readline().split(':')[1])
    for _ in range(3):
      line = streamIn.readline()
    self.calcSetup['numberOfFrozenCoreElectrons'] = int(streamIn.readline().split(':')[1])
    line = streamIn.readline()
    self.calcSetup['numberOfFrozenVirtualOrbitals'] = int(streamIn.readline().split(':')[1])
    self.calcSetup['molecularSpinMultiplicity'] = int(streamIn.readline().split(':')[1])
    while line.find('Correlation Information')<0:
      line = streamIn.readline()
    line = streamIn.readline()
    self.calcSetup['waveFunctionTheory'] = streamIn.readline().split(':')[1]
    perturbative = self.calcSetup['waveFunctionTheory'].find('perturbation')>=0
    if perturbative:
      self.calcSetup['waveFunctionTheory'].rstrip('w/ perturbation')
      self.calcSetup['waveFunctionTheory']+=' with '+streamIn.readline().split(':')[1]+' perturbative correction'
    self.calcTask['calculationSetup'] = self.calcSetup
    while line:
      if line.find('correlation energy')>=0:
        vars = line.split('=')
        theory = vars[0].split('correlation')[0]
        self.calcRes['correlationEnergy'].append({'theory': theory, 'value' : float(vars[1]), 'units' : 'Hartree'})
      elif line.find('correction energy')>=0:
        vars = line.split('=')
        theory = vars[0].split('correlation')[0]
        self.calcRes['perturbativeCorrectionEnergy'].append({'theory': theory, 'value' : float(vars[1]), 'units' : 'Hartree'})
      elif line.find('total energy')>=0:
        vars = line.split('=')
        theory = vars[0].split('correlation')[0]
        self.calcRes['totalEnergy'].append({'theory': theory, 'value' : float(vars[1]), 'units' : 'Hartree'})
      elif line.find('Task  times')>=0:
        self.calcTask['calculationTime'] = self.readTaskTimes(line)
        break
      line = streamIn.readline()
    self.calcTask['calculationResults'] = self.calcRes
    self.calculations.append(self.calcTask)
  
class basisObj:
  def __init__(self):
    self.basCount = 0
    self.basUpdated = False
    self.basis = {}
  def readBasis(self,line,streamIn):
    self.basCount += 1
    self.basis = {}
    atomBasCount = 0
    if not self.basUpdated : 
      self.basis = {}
      self.basis['id'] = 'BasisSet.'+str(self.basCount)
      self.basis['basisFunctions'] = []
    cartSpher = line.replace('(','').replace(')','').split()[-1]
    for _ in range(2): 
      line = streamIn.readline()
    while line:
      if line.find('Summary of "ao')>=0:
        count = 0
        for _  in range(4): 
          line = streamIn.readline()
        while line:
          vars = line.split()
          if len(vars) == 0:
            break
          else:
            self.basis['basisFunctions'][count].update({'basisSetName': vars[1]})
            count += 1
          line = streamIn.readline()
        break
      atomLab, atomName = line.replace('(','').replace(')','').split()
      atomBasCount += 1
      basisAtom = {}
      basisAtom['id'] = atomLab+"-orb"+str(atomBasCount)
      basisAtom['elementLabel'] = atomLab
      basisAtom['elementName'] = atomName
      basisAtom['basisSetType'] = 'orbitalBasis'
      basisAtom['basisSetHarmonicType'] = cartSpher
      for _ in range(4): 
        line = streamIn.readline()
      cont = 1
      basisExp = []
      basisCoef = []
      basisCont = {}
      basisConts = []
      while line:
        vars = line.split()
        if len(vars) != 0 and len(vars) != 4:
          basisCont['basisSetExponent'] = basisExp
          basisCont['basisSetCoefficient'] = basisCoef
          basisConts.append(basisCont)
          break
        if len(vars) == 4:
          if int(vars[0]) == cont:
            basisCont['id'] = atomLab+"-orbc"+str(cont)
            basisCont['basisSetShellType'] = vars[1].lower()
            basisExp.append(vars[2])
            basisCoef.append(vars[3])
          else:
            basisCont['basisSetExponent'] = basisExp
            basisCont['basisSetCoefficient'] = basisCoef
            basisConts.append(basisCont)
            cont+=1
            basisExp = []
            basisCoef = [] 
            basisCont = {}
            basisCont['id'] = atomLab+"-orbc"+str(cont)
            basisCont['basisSetShellType'] = vars[1].lower()
            basisExp.append(float(vars[2]))
            basisCoef.append(float(vars[3]))
        line = streamIn.readline()
      basisAtom['basisSetContraction'] = basisConts
      self.basis['basisFunctions'].append(basisAtom)
    self.basUpdated = True
  def readEcp(self,line,streamIn):
    atomBasCount = 0
    cartSpher = line.replace('(','').replace(')','').split()[-1]
    if not self.basUpdated : 
      self.basis = {}
      self.basis['id'] = 'BasisSet.'+str(self.basCount)
      self.basis['basisFunctions'] = []
    line = streamIn.readline()
    emptyLine = 0
    while line:
      if line.find('Module')>=0 or line.find('library')>=0 or line.find('NWChem')>=0 or emptyLine>1: 
        break
      vars = streamIn.readline().split()
      atomLab, atomName, elec = vars[0], vars[1].replace('(','').replace(')',''), vars[3]
      atomBasCount += 1
      ecpAtom = {}
      ecpAtom['id'] = atomLab+"-ecp"+str(atomBasCount)
      ecpAtom['elementLabel'] = atomLab
      ecpAtom['elementName'] = atomName
      ecpAtom['basisSetType'] = 'ecpBasis'
      ecpAtom['numberElectronsReplaced'] = elec
      ecpAtom['basisSetHarmonicType'] = cartSpher
      ecpExp = []
      ecpCoef = [] 
      ecpRExp = [] 
      ecpCont = {}
      ecpConts = []
      for _ in range(4): 
        line = streamIn.readline()
      cont = 1
      while line:
        if line.find('Module')>=0 or line.find('library')>=0 or line.find('electrons')>=0: 
          break
        vars = line.split()
        if len(vars) == 0:
          emptyLine += 1
        if (len(vars) != 0 and len(vars) != 5) or emptyLine>1:
          break
        if len(vars) == 5:
           emptyLine = 0
           if int(vars[0]) == cont:
             ecpCont['id'] = atomLab+"-ecpc"+str(cont)
             ecpCont['basisSetShellType'] = vars[1].upper()
             ecpRExp.append(vars[2])
             ecpExp.append(vars[3])
             ecpCoef.append(vars[4])
           else:
             ecpCont['basisSetRExponent'] = ecpRExp
             ecpCont['basisSetExponent'] = ecpExp
             ecpCont['basisSetCoefficient'] = ecpCoef
             ecpConts.append(ecpCont)
             cont+=1
             ecpExp = []
             ecpCoef = [] 
             ecpRExp = [] 
             ecpCont = {}
             ecpCont['id'] = atomLab+"-ecpc"+str(cont)
             ecpCont['basisSetShellType'] = vars[1].upper()
             ecpRExp.append(float(vars[2]))
             ecpExp.append(float(vars[3]))
             ecpCoef.append(float(vars[4]))
        line = streamIn.readline()
      ecpAtom['basisSetContraction'] = ecpConts
      self.basis['basisFunctions'].append(ecpAtom)
    self.basUpdated = True
    
class moleculeObj:
  def __init__(self):
    self.molCount = 0
    self.geomUpdated = False
    self.molecule = {}
    self.atomCount = 0
    self.molecularFormula = ''
  def readGeom(self,line,streamIn):
    elements = ['blank','H','He',
                'Li','Be','B','C','N','O','F','Ne',
                'Na','Mg','Al','Si','P','S','Cl','Ar',
                'K','Ca','Sc','Ti','V','Cr','Mn','Fe',
                'Co','Ni','Cu','Zn','Ga','Ge','As',
                'Se','Br','Kr','Rb','Sr','Y','Zr','Nb',
                'Mo','Tc','Ru','Rh','Pd','Ag','Cd',
                'In','Sn','Sb','Te','I','Xe','Cs','Ba',
                'La','Ce','Pr','Nd','Pm','Sm','Eu',
                'Gd','Tb','Dy','Ho','Er','Tm','Yb',
                'Lu','Hf','Ta','W','Re','Os','Ir','Pt',
                'Au','Hg','Tl','Pb','Bi','Po','At',
                'Rn','Fr','Ra','Ac','Th','Pa','U','Np',
                'Pu','Am','Cm','Bk','Cf','Es','Fm',
                'Md','No','Lr','Rf','Db','Sg','Bh',
                'Hs','Mt','Ds','Rg','Cn','Uut','Fl',
                'Uup','Lv','Uus','Uuo']
    self.molCount += 1
    self.molecule = {}
    self.atomCount = 0
    self.molecule['id'] = 'Molecule.'+str(self.molCount)
    self.molecularFormula = ''
    atoms = []
    for _ in range(3): 
      line = streamIn.readline()
    if line.split()[3] == 'a.u.':
      geomUnit = 'bohr'
    else:
      geomUnit = 'angstrom'
    for _ in range(4): 
      line = streamIn.readline()
    formulaList=[]
    while line:
      vars = line.split()
      if len(vars) < 6:
        break
      else:
        atom = {}
        cart = {} 
        val = []
        atom['id'] = 'Atom.'+str(vars[0])+'.Mol.'+str(self.molCount)
        atom['elementLabel'] = vars[1]
        elementNumber = int(float(vars[2]))
        atom['elementNumber'] = elementNumber
        atom['elementSymbol'] = elements[elementNumber]
        if any(elements[elementNumber] in element for element in formulaList):
          for element in formulaList:
            if elements[elementNumber] in element:
              element[1] += 1
        else:
          formulaList.append([elements[elementNumber],1])
        val = [float(vars[3]),float(vars[4]),float(vars[5])]
        cart['value'] = val
        cart['units'] = geomUnit
        atom['cartesianCoordinates'] = cart
        atoms.append(atom)
        self.atomCount += 1
      line = streamIn.readline()
    for item in formulaList:
      if item[0] == 'C':
        if item[1] == 1:
          self.molecularFormula += item[0]
        else:
          self.molecularFormula += item[0]+str(item[1])
        formulaList.remove(item)
        break
    for item in formulaList:
      if item[0] == 'H':
        if item[1] == 1:
          self.molecularFormula += item[0]
        else:
          self.molecularFormula += item[0]+str(item[1])
        formulaList.remove(item)
        break
    formulaList.sort(key=lambda x: x[0])
    for item in formulaList:
      if item[1] == 1:
        self.molecularFormula += item[0]
      else:
        self.molecularFormula += item[0]+str(item[1])
    self.molecule['atoms'] = atoms
    symmetry = {}
    symmetry['groupname'] = 'C1'
    while line:
      if line.find('Module')>=0 or line.find('library'): 
        break
      if line.find('Symmetry info')>=0:
        for _ in range(3): 
          line = streamIn.readline()
        symmetry['groupName'] = line.split()[2]
        break
      line = streamIn.readline()
    self.molecule['symmetry'] = symmetry
    self.geomUpdated = True
