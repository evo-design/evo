from setuptools import setup, find_packages

with open('evo/version.py') as infile:
    exec(infile.read())

with open('README.md') as f:
    readme = f.read()

with open('requirements.txt') as f:
    requirements = f.read().split('\n')

sources = {
    'evo': 'evo',
    'evo.scripts': 'scripts',
}

setup(
    name='evo-model',
    version=version,
    description='DNA foundation modeling from molecular to genome scale.',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Team Evo',
    url='https://github.com/evo-design/evo',
    license='Apache-2.0',
    packages=sources.keys(),
    package_data={'evo': ['evo/configs/*.yml']},
    include_package_data=True,
    package_dir=sources,
    install_requires=requirements,
    python_requires='>=3.6',
)
