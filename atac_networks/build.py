from setuptools.command.build_ext import build_ext
from setuptools import Extension
import numpy
import platform


if platform.system() == "Darwin":
    extra_compile_args = ["-I/System/Library/Frameworks/vecLib.framework/Headers"]
    if "ppc" in platform.machine():
        extra_compile_args.append("-faltivec")

    extra_link_args = ["-Wl,-framework", "-Wl,Accelerate"]
    include_dirs = [numpy.get_include()]

else:
    include_dirs = [numpy.get_include(), "/usr/local/include"]
    extra_compile_args = ["-msse2", "-O2", "-fPIC", "-w"]
    extra_link_args = ["-llapack"]

from distutils.errors import DistutilsPlatformError, CCompilerError, DistutilsExecError, DistutilsPlatformError


ext_modules = [
    Extension(name="atac_networks.pyquic",
              include_dirs=[numpy.get_include()],
              sources=["atac_networks/pyquic/QUIC.C", "atac_networks/pyquic/pyquic.pyx"],
              extra_compile_args=extra_compile_args,
              extra_link_args=extra_link_args,
              language="c++",
             ),
]


class BuildFailed(Exception):
    pass


class ExtBuilder(build_ext):
    def run(self):
        try:
            build_ext.run(self)
        except (DistutilsPlatformError, FileNotFoundError):
            pass

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsExecError, DistutilsPlatformError, ValueError):
            pass


def build(setup_kwargs):
    """
    This function is mandatory in order to build the extensions.
    """
    setup_kwargs.update(
        {"ext_modules": ext_modules, "cmdclass": {"build_ext": ExtBuilder}}
    )