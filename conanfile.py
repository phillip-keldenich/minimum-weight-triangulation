from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.files import load, copy
import os, subprocess


class CGALMWTConan(ConanFile):
    name = "cgal_mwt"

    def set_version(self):
        """
        Set the version from the file 'version.txt'.
        If the command line gives another value by
        the flag --version, that value is used instead.
        """
        self.version = load(self, os.path.join(self.recipe_folder, "version.txt")).strip()

    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False], "testing": [True, False]}
    default_options = {"shared": False, "fPIC": True, "testing": False}

    def config_options(self):
        """
        Routine checks which options are available,
        e.g., depending on OS, version and compiler.
        """
        if self.settings.get_safe("os") == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        """
        Routine that checks options for validity.
        If, e.g., shared builds are unavailable on Windows,
        we can raise a ConanInvalidConfiguration error here.
        """
        pass

    def build_requirements(self):
        """
        Method that declares our build requirements.
        """
        try:
            # try to avoid pulling in an extra copy of cmake if the system already has one
            subprocess.run(["cmake", "--help"], text=True, capture_output=True, check=True)
        except Exception:
            self.tool_requires("cmake/[>=3.16]")

    def requirements(self):
        self.requires("cgal/[>=5.6]")
        self.requires("onetbb/[>=2021.3.0]")
        self.requires("nlohmann_json/[>=3.9.1]")
        self.requires("gurobi/[>=10]@ibralg/develop")
        self.requires("boost/[>=1.76.0]")
        if self.options.testing:
            self.requires("doctest/[>=2.4.3]")

    def export(self):
        """
        Copy anything the conanfile.py needs to run
        to self.export_folder. Also, if there is a 
        special license for just the conanfile.py,
        that should be copied here as well.
        """
        copy(self, "version.txt", self.recipe_folder, self.export_folder)

    def export_sources(self):
        """
        Copy sources (useful if the conanfile.py is in the same repo as the sources).
        """
        copy(self, "src/*", self.recipe_folder, self.export_sources_folder, keep_path=True)
        copy(self, "include/CGAL_MWT/*", self.recipe_folder, self.export_sources_folder, keep_path=True)
        copy(self, "CMakeLists.txt", self.recipe_folder, self.export_sources_folder, keep_path=True)
        copy(self, "*.md", self.recipe_folder, self.export_sources_folder, keep_path=True)
        copy(self, "version.txt", self.recipe_folder, self.export_sources_folder, keep_path=True)

    def source(self):
        """
        Could download sources (or additional sources) here.
        Useful if the sources are in a Git other than 
        the git hosting the conanfile.py.
        """
        pass

    def layout(self):
        """
        There is usually no need to change this.
        """
        cmake_layout(self)

    def generate(self):
        """
        Call CMake, possibly setting CMake variables.
        """
        tc = CMakeToolchain(self)
        # dictionary of CMake variables to set
        cmake_vars_to_set = {
            "CGAL_MWT_ENABLE_TESTING": "ON" if self.options.testing else "OFF"
        }
        for var, val in cmake_vars_to_set.items():
            tc.variables[var] = val
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        """
        Call cmake --build.
        """
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        """
        Call cmake --install.
        """
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        """
        Describe libraries we export, their inter-dependencies,
        their dependencies on system libraries and other libraries,
        additional flags needed by consumers, ...
        """
        to_libname = lambda name: f"{name}.lib" if self.settings.os == "Windows" else f"lib{name}.a"

        self.cpp_info.components["cgal_mwt"].libs = [to_libname("cgal_mwt")]
        self.cpp_info.components["cgal_mwt"].requires = ["cgal::cgal", "onetbb::onetbb", 
                                                         "nlohmann_json::nlohmann_json",
                                                         "gurobi::gurobi",
                                                         "boost::boost"]
        if self.settings.compiler == "gcc":
            self.cpp_info.components["cgal_mwt"].cflags = ["-frounding-math"]
            self.cpp_info.components["cgal_mwt"].cxxflags = ["-frounding-math"]
        elif self.settings.compiler == "msvc":
            self.cpp_info.components["cgal_mwt"].cflags = ["/fp:strict", "/fp:except-"]
            self.cpp_info.components["cgal_mwt"].cxxflags = ["/fp:strict", "/fp:except-"]
