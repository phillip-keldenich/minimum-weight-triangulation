import subprocess
try:
    # conan v2
    from conan import ConanFile
    from conan.tools.cmake import CMake, CMakeToolchain, cmake_layout    
    is_conan_v2 = True
except ImportError:
    # conan v1
    from conans import ConanFile, tools
    from conan.tools.cmake import CMake, CMakeToolchain
    is_conan_v2 = False


class GurobiTestPackageConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps", "VirtualBuildEnv", "VirtualRunEnv"
    options = {}

    def build_requirements(self):
        """
        Method that declares our build requirements.
        """
        try:
            # try to avoid pulling in an extra copy of cmake if the system already has one
            subprocess.run(["cmake", "--help"], text=True, capture_output=True, check=True)
        except Exception:
            self.tool_requires("cmake/[>=3.16]")

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def test(self):
        cmake = self._configure_cmake()
        cmake.test()


def _requires_method(self):
    self.requires(self.tested_reference_str)

def _layout_method(self):
    cmake_layout(self)


if is_conan_v2:
    GurobiTestPackageConan.requirements = _requires_method
    GurobiTestPackageConan.layout = _layout_method
