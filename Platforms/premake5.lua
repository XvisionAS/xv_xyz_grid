local action = _ACTION or ""

workspace "xv_xyz_grid"
	location 		("../build/" .. action)
	targetdir		("../build/bin")

  configurations {"Debug", "Release"}
	defines { "GLM_FORCE_RADIANS" }
	defines     { "TRILIBRARY", "NO_TIMER", "TRILIBRARY_x64", "ANSI_DECLARATORS" }
	filter "action:vs*"
		defines     { "_CRT_SECURE_NO_WARNINGS" }

  project "xv_xyz_grid"
    kind "ConsoleApp"


		language "C++"
		files {
			"../Externals/triangle/**.cc",
			"../Externals/triangle/**.h",
			"../Externals/glm/**.hpp",
			"../Externals/cmdline/**.h",			
			"../*.h",
			"../*.cpp"
		}

		includedirs {
			"../Externals",
		}

		filter { "files:../*.cpp" }

		flags { "C++11" }
		filter "system:linux"
			links {"m", "dl"}
    	buildoptions { "-std=c++11", "-m64" }
			linkoptions { "-m64" }

		filter "system:macosx"
			platforms 	{"Universal"}
			defines {"MACOSX"}
			links       { "CoreFoundation.framework", "Cocoa.framework", "IOKit.framework"}
			buildoptions { "-mmacosx-version-min=10.4" }
			linkoptions  { "-mmacosx-version-min=10.4" }

		filter "system:windows"
			architecture ("x86_64")			

		filter "configurations:Debug"
			defines { "DEBUG" }
			flags { "Symbols" }
		filter "configurations:Release"
			defines { "NDEBUG" }
			optimize "Full"
