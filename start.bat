@echo off
REM Drosophila Research Assistant - Windows Startup Script

echo ================================================
echo  Drosophila Research Assistant
echo ================================================
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python is not installed or not in PATH
    echo Please install Python 3.8 or higher from python.org
    pause
    exit /b 1
)

echo [OK] Python found
echo.

REM Check if virtual environment exists
if not exist "venv\" (
    echo [SETUP] Creating virtual environment...
    python -m venv venv
    if errorlevel 1 (
        echo [ERROR] Failed to create virtual environment
        pause
        exit /b 1
    )
    echo [OK] Virtual environment created
)

REM Activate virtual environment
echo [SETUP] Activating virtual environment...
call venv\Scripts\activate.bat

REM Check if requirements are installed
if not exist "venv\.installed" (
    echo.
    echo [SETUP] Installing dependencies (this may take a minute)...
    python -m pip install --upgrade pip
    pip install -r requirements.txt
    if errorlevel 1 (
        echo [ERROR] Failed to install dependencies
        pause
        exit /b 1
    )
    echo. > venv\.installed
    echo [OK] Dependencies installed
)

REM Check for API key
if "%ANTHROPIC_API_KEY%"=="" (
    if exist ".env" (
        echo [OK] Loading API key from .env file...
        for /f "tokens=1,2 delims==" %%a in (.env) do set %%a=%%b
    ) else (
        echo.
        echo [WARNING] ANTHROPIC_API_KEY not found!
        echo.
        echo Please set your API key using one of these methods:
        echo.
        echo Option 1 - Set for this session:
        echo   set ANTHROPIC_API_KEY=your_api_key_here
        echo   start.bat
        echo.
        echo Option 2 - Create a .env file:
        echo   echo ANTHROPIC_API_KEY=your_api_key_here ^> .env
        echo   start.bat
        echo.
        echo Get your API key from: https://console.anthropic.com
        echo.
        set /p api_key="Enter your Anthropic API key (or press Enter to exit): "
        if "%api_key%"=="" (
            echo Exiting. Please set your API key and try again.
            pause
            exit /b 1
        )
        set ANTHROPIC_API_KEY=%api_key%
        echo ANTHROPIC_API_KEY=%api_key% > .env
        echo [OK] API key saved to .env file
    )
)

echo [OK] API key configured
echo.
echo ================================================
echo  Starting Drosophila Assistant
echo ================================================
echo.
echo  Open in your browser: http://localhost:5000
echo.
echo  Press Ctrl+C to stop the server
echo.
echo ================================================
echo.

REM Start the application
python app.py

pause
