from distutils.core import setup

setup(
      name='DpModel',
      version='1.0.0',
      package_dir = {'': 'lib'},
      packages=[
                'dp_model',
                'dp_model.post_process',
                'dp_model.post_process.plot'
                ],
      scripts=['DpModel'],
     )
