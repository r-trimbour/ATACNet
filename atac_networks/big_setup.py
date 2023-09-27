# -*- coding: utf-8 -*-
from setuptools import setup

packages = \
['atac_networks', 'atac_networks.pyquic']

package_data = \
{'': ['*']}

install_requires = \
['joblib>=1.3.2,<2.0.0',
 'numpy>=1.26.0,<2.0.0',
 'pandas>=2.1.1,<3.0.0',
 'scikit-learn>=1.3.1,<2.0.0',
 'scipy>=1.11.2,<2.0.0',
 'setuptools>=68.2.2,<69.0.0',
 'tqdm>=4.66.1,<5.0.0']

setup_kwargs = {
    'name': 'atac-networks',
    'version': '0.1.0',
    'description': '',
    'long_description': '',
    'author': 'r-trimbour',
    'author_email': 'trimbour@edu.bio.ens.psl.eu',
    'maintainer': 'None',
    'maintainer_email': 'None',
    'url': 'None',
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'python_requires': '>=3.9,<3.13',
}
from build import *
build(setup_kwargs)

setup(**setup_kwargs)

