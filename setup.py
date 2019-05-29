import setuptools


def readme():
    with open('README.rst') as f:
        return f.read()

setuptools.setup(name='word-algebras-sage',
      #python_requires='>=3.5.2',
      version='0.1.0',
      packages=['word_algebras'],
      description='A SAGE library implementing the (nilpotent) shuffle and concatenation algebras.',
      long_description=readme(),
      author='Joscha Diehl',
      #author_email='',
      url='https://github.com/diehlj/word-algebras-sage',
      license='Eclipse Public License',
      install_requires=['numpy', 'scipy', 'sympy'],
      #setup_requires=['setuptools_git >= 0.3', ],
      #test_suite='tests'
      )
