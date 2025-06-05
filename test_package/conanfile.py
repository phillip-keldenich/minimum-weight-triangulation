from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout
from conan.tools.build import can_run


class CGALMWTTestPackageConan(ConanFile):
      settings = "os", "compiler", "build_type", "arch"
      generators = "CMakeDeps", "CMakeToolchain"

      def requirements(self):
            self.requires(self.tested_reference_str)
            self.requires("doctest/[>=2.4.3]")
            self.requires("cgal/[>=5.6]")
            self.requires("onetbb/[>=2021.3.0]")
            self.requires("nlohmann_json/[>=3.9.1]")
            self.requires("gurobi/[>=10]@ibralg/develop")
            self.requires("boost/[>=1.76.0]")

      def build(self):
            cmake = CMake(self)
            cmake.configure()
            cmake.build()

      def layout(self):
            cmake_layout(self)

      def test(self):
            if can_run(self):
                  cmake = CMake(self) 
                  cmake.test()
