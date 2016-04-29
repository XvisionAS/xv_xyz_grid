local action = _ACTION or ""

workspace "xv_xyz_grid"
	location 		("../build/" .. action)
	targetdir		("../build/bin")

  configurations {"Debug", "Release"}
	defines { "GLM_FORCE_RADIANS" }
	defines     { "SINGLE", "TRILIBRARY", "NO_TIMER" }
	filter "action:vs*"
		defines     { "_CRT_SECURE_NO_WARNINGS" }

  project "xv_xyz_grid"
    kind "ConsoleApp"
    language "C"
		files {
			"../Externals/triangle/**.c",
			"../Externals/triangle/**.h"
		}
		language "C++"
		files {
			"../Externals/glm/**.hpp",
			"../Externals/cmdline/**.h",
			"../xv_xyz_grid.cpp"
		}
		includedirs {
			"../Externals",
		}
		filter "system:linux"
			toolset "clang"
			links {"m", "dl"}
		filter "system:macosx"
			platforms 	{"Universal"}
			defines {"MACOSX"}
			links       { "CoreFoundation.framework", "Cocoa.framework", "IOKit.framework"}
			buildoptions { "-mmacosx-version-min=10.4" }
			linkoptions  { "-mmacosx-version-min=10.4" }
		filter "system:windows"
			architecture ("x86")
		filter "configurations:Debug"
			defines { "DEBUG" }
			flags { "Symbols" }
		filter "configurations:Release"
			defines { "NDEBUG" }
			optimize "On"
