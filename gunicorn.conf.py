# Gunicorn configuration file for Drosophila Research Assistant
# Place this file in your project root as: gunicorn.conf.py

import multiprocessing
import os

# Server socket
bind = f"0.0.0.0:{os.environ.get('PORT', '5000')}"
backlog = 2048

# Worker processes
workers = 2  # Don't use too many on free tier
worker_class = 'sync'
worker_connections = 1000
timeout = 120  # 2 minutes - crucial for Claude API calls
keepalive = 5
max_requests = 1000  # Restart workers after 1000 requests to prevent memory leaks
max_requests_jitter = 50

# Logging
accesslog = '-'  # Log to stdout
errorlog = '-'   # Log to stderr
loglevel = 'info'
access_log_format = '%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s" %(D)s'

# Process naming
proc_name = 'drosophila_assistant'

# Server mechanics
daemon = False
pidfile = None
umask = 0
user = None
group = None
tmp_upload_dir = None

# Performance tuning
preload_app = False  # Don't preload - can cause issues with API keys
worker_tmp_dir = '/dev/shm'  # Use RAM disk for worker temp files (faster)

# Graceful timeout
graceful_timeout = 30

print("="*60)
print("🪰 Gunicorn Configuration Loaded")
print("="*60)
print(f"Workers: {workers}")
print(f"Timeout: {timeout}s")
print(f"Port: {os.environ.get('PORT', '5000')}")
print("="*60)