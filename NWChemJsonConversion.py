from nwchemToJson import *
for files in range(1,len(sys.argv)):
  fileIn = open(sys.argv[files],'r')
  fileOut = open(sys.argv[files]+'.json','w')
  jsonObj = nwchemToJson('noOrbitals')
  fileOut.write(jsonObj.convert(fileIn))
  fileIn.close()
  fileOut.close()
