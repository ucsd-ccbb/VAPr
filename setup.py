import sys
from setuptools import setup


def setup_package():
    needs_sphinx = {'build_sphinx', 'upload_docs'}.intersection(sys.argv)
    sphinx = ['sphinx'] if needs_sphinx else []
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
          setup_requires=['six'] + sphinx)


if __name__ == "__main__":
    setup_package()
