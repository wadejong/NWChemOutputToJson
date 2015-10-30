from nwchemToJson import *
if sys.argv[1] == 'noOrbitals':
  start = 2
  argument = 'noOrbitals'
else:
  start = 1
  argument = ''
for files in range(start,len(sys.argv)):
  fileIn = open(sys.argv[files],'r')
  fileOut = open(sys.argv[files]+'.json','w')
  jsonObj = nwchemToJson(argument)
  fileOut.write(jsonObj.convert(fileIn))
  fileIn.close()
  fileOut.close()
