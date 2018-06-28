##############################################################################
# This source file is part of the NWChemOutputToJson project.
# Copyright (c) 2016, The Regents of the University of California, through 
# Lawrence Berkeley National Laboratory (subject to receipt of any required 
# approvals from the U.S. Dept. of Energy).
# This source code is released under the BSD 3-Clause License, (the "License").
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

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
  print('Converting file ',sys.argv[files])
  jsonObj = nwchemToJson(argument)
  fileOut.write(jsonObj.convert(fileIn))
  fileIn.close()
  fileOut.close()
