@REM @echo off

@REM echo Building and running project using Visual Studio in release mode
@REM cd build && cmake ..\ && MSBuild.exe src\friction_dynamics.vcxproj /p:Configuration=Release && src\Release\friction_dynamics.exe

@echo off

if "%1" == "" (
    echo Please specify compiler [msvc, mingw]
    exit
)

if "%1" == "msvc" (
    if not exist build (
        mkdir build
    )

    if "%2" == "tests" (
        echo Building and running all tests using Visual Studio
        cd build && cmake .. -G "Visual Studio 19 2017 Win64" && MSBuild.exe tests.vcxproj && cd Debug && tests.exe
        exit
        
    ) else (
        if "%2" == "release" (
            echo Building and running project using Visual Studio in release mode
            cd build && cmake .. -DCMAKE_GENERATOR_PLATFORM=x64 && MSBuild.exe src\friction_dynamics.vcxproj /p:Configuration=Release && cd src\Release && friction_dynamics.exe
            exit
        ) else (
            echo Building and running project using Visual Studio
            cd build && cmake .. && MSBuild.exe src\friction_dynamics.vcxproj && cd src\Debug && friction_dynamics.exe
            exit
        )
    )
)

if "%1" == "mingw" (
    if not exist build\mingw (
        mkdir build\mingw
    )

    if "%2" == "tests" (
        echo Building and running all tests using MinGW Makefiles (SLOW with Catch2!!) && cd build\mingw && cmake ..\..\ -G "MinGW Makefiles" && make tests -j2 && .\tests.exe
    ) else "%2" == "" (
        echo Building and running project using MinGW Makefiles &&  cd build\mingw && cmake ..\..\ -G "MinGW Makefiles" && make critplane -j2 && .\fr.exe
    )
) 