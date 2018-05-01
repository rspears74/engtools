from setuptools import setup

setup(name='engtools',
      version='0.1',
      description='Structural engineering related functions and tools.',
      url='https://github.com/rspears74/engtools',
      author='Randall Spears',
      author_email='rspears690@gmail.com',
      license='MIT',
      packages=['engtools'],
      install_requires=[
          'prettytable',
      ],
      zip_safe=False)