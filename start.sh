#!/bin/bash

# Fly Design - Startup Script
# This script makes it easy to start the application

echo "🪰 Fly Design"
echo "================================"
echo ""

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.11."
    exit 1
fi

echo "✓ Python found: $(python3 --version)"

# Check if virtual environment exists
if [ ! -d "venv" ]; then
    echo ""
    echo "📦 Creating virtual environment..."
    python3 -m venv venv
    
    if [ $? -ne 0 ]; then
        echo "❌ Failed to create virtual environment"
        exit 1
    fi
    
    echo "✓ Virtual environment created"
fi

# Activate virtual environment
echo ""
echo "🔧 Activating virtual environment..."
source venv/bin/activate

# Check if requirements are installed
if [ ! -f "venv/.installed" ]; then
    echo ""
    echo "📥 Installing dependencies (this may take a minute)..."
    pip install --upgrade pip
    pip install -r requirements.txt
    
    if [ $? -ne 0 ]; then
        echo "❌ Failed to install dependencies"
        exit 1
    fi
    
    touch venv/.installed
    echo "✓ Dependencies installed"
fi

# Check for API key
if [ -z "$ANTHROPIC_API_KEY" ] && [ ! -f ".env" ]; then
    echo ""
    echo "⚠️  ANTHROPIC_API_KEY not found!"
    echo ""
    echo "Please set your API key using one of these methods:"
    echo ""
    echo "Option 1 - Set for this session:"
    echo "  export ANTHROPIC_API_KEY=your_api_key_here"
    echo "  ./start.sh"
    echo ""
    echo "Option 2 - Create a .env file:"
    echo "  echo 'ANTHROPIC_API_KEY=your_api_key_here' > .env"
    echo "  ./start.sh"
    echo ""
    echo "Get your API key from: https://console.anthropic.com"
    echo ""
    read -p "Do you want to enter your API key now? (y/n) " -n 1 -r
    echo ""
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        read -p "Enter your Anthropic API key: " api_key
        export ANTHROPIC_API_KEY=$api_key
        echo "ANTHROPIC_API_KEY=$api_key" > .env
        echo "✓ API key saved to .env file"
    else
        echo "Exiting. Please set your API key and try again."
        exit 1
    fi
fi

echo ""
echo "✓ Configuration check complete"
echo ""
echo "================================"
echo "🚀 Starting Fly Design"
echo "================================"
echo ""
echo "🌐 Open in your browser: http://localhost:5000"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Start the application
python3 app.py
