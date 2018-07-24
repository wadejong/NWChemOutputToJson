##############################################################################
# This source file is part of the NWChemOutputToJson project.
# Copyright (c) 2018, The Regents of the University of California, through 
# Lawrence Berkeley National Laboratory (subject to receipt of any required 
# approvals from the U.S. Dept. of Energy).
# This source code is released under the BSD 3-Clause License, (the "License").
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

from distutils.core import setup
import setuptools

setup(
    name='nwchem2json',
    version='1.0.0',
    url='https://github.com/wadejong/NWChemOutputToJson',
    author='Bert de Jong',
    description='Python files for reading NWChem output and converting to Json',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'nwchem2json=nwchem2json.scripts.nwchem2json:main',
        ],
    },
)
