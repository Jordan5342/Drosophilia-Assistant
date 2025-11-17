#!/usr/bin/env python3
"""Test the optimized Drosophila Assistant with FlyBase publications"""

import os
from drosophila_assistant.py import DrosophilaAssistant

def test_assistant():
    """Test that FlyBase publications are being fetched"""
    
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        print("❌ ANTHROPIC_API_KEY not set")
        return
    
    assistant = DrosophilaAssistant(api_key)
    
    print("\n" + "="*70)
    print("TEST 1: InR Gene (should find 13+ PMIDs from FlyBase)")
    print("="*70)
    response = assistant.chat("Tell me about the InR gene")
    print("\nAssistant response (first 500 chars):")
    print(response[:500] + "...\n")
    
    print("\n" + "="*70)
    print("TEST 2: Notch Gene (should find 23+ PMIDs from FlyBase)")
    print("="*70)
    assistant.reset_conversation()
    response = assistant.chat("What is the Notch gene in Drosophila?")
    print("\nAssistant response (first 500 chars):")
    print(response[:500] + "...\n")
    
    print("\n" + "="*70)
    print("TEST 3: foxo Gene (0 PMIDs on FlyBase, should use PubMed fallback)")
    print("="*70)
    assistant.reset_conversation()
    response = assistant.chat("Tell me about foxo")
    print("\nAssistant response (first 500 chars):")
    print(response[:500] + "...\n")

if __name__ == "__main__":
    test_assistant()