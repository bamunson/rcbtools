from setuptools import setup

setup(name='rcbtools',
      version='0.1',
      description='Analysis toolkit for RCB stars evolved in MESA',
      url='https://github.com/bamunson/rcbtools',
      author='Brad Munson',
      author_email='brad.munson94@gmail.com',
      license='MIT',
      packages=['rcbtools'],
      install_requires=['numpy',
			'matplotlib',
			'pandas',                     
                       ],

)