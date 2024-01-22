from distutils.core import setup
setup(
  name = 'MBN_tools',
  packages = ['MBN_tools'],
  version = '0.1',
  license='MIT',
  description = 'Useful tools for MBN Explorer simulations and analysis',
  author = 'Matthew Dickers',
  author_email = 'mattdickers@gmail.com',
  url = 'https://github.com/mattdickers/MBN_tools',
  download_url = 'https://github.com/mattdickers/MBN_tools/archive/refs/tags/v_0.1.tar.gz',
  keywords = ['MBN Explorer', 'MBN Studio', 'Molecular Dynamics', 'MD', 'Data Analysis'],
  install_requires=[
          'numpy',
          'mdtraj',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
  ],
)