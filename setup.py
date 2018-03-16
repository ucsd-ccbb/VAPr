from setuptools import setup

setup(name='VAPr',
      version='3.0',
      description='Package for NoSQL variant data storage, annotation and prioritization.',
      url='https://github.com/ucsd-ccbb/VAPr',
      author='Carlo Mazzaferro, Amanda Birmingham, Adam Mark',
      author_email='cmazzafe@ucsd.edu',
      install_requires=['pymongo', 'myvariant', 'pyvcf', 'pandas', 'tqdm'],
      license='MIT',
      packages=['VAPr'],
      zip_safe=False,
      extras_require={
            'docs': [
                  'sphinx >= 1.7',
                  'sphinx_rtd_theme']})
