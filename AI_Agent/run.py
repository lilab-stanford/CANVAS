#!/usr/bin/env python3
"""
Quick runner for Spatial Feature Agent
"""

from spatial_agent import SpatialFeatureAgent
import os
import getpass

def quick_demo():
    """Run a quick demonstration of the agent"""
    # Header
    print("╔" + "═" * 70 + "╗")
    print("║" + " " * 14 + "AI Agent for Habitat-Level Interpretation" + " " * 15 + "║")
    print("║" + " " * 22 + "of CANVAS Spatial Features" + " " * 22 + "║")
    print("╚" + "═" * 70 + "╝")
    print()
    
    # System status
    print("📡 SYSTEM STATUS")
    print("─" * 50)
    api_key = os.getenv("OPENAI_API_KEY")
    if api_key:
        print("🔑 API Configuration: ✓ Found OPENAI_API_KEY")
    else:
        print("🔑 API Configuration: ✗ OPENAI_API_KEY not set")
        try:
            api_key = getpass.getpass("Enter OpenAI API key: ").strip()
        except (KeyboardInterrupt, EOFError):
            api_key = ""
        if not api_key:
            print("⚠️ No API key provided. Exiting quick demo.")
            return
        print("🔑 API Configuration: ✓ Received from prompt")

    api_base_url = os.getenv("OPENAI_BASE_URL") or os.getenv("OPENAI_API_BASE")
    if api_base_url:
        print(f"🌐 API Endpoint Override: {api_base_url}")
    else:
        print("🌐 API Endpoint: OpenAI default (override via OPENAI_BASE_URL if needed)")

    print("🔧 Agent Initialization: ", end="")

    try:
        agent = SpatialFeatureAgent(
            api_key=api_key,
            model_name="gpt-4o",
            temperature=0.3,
            num_workers=4,
            api_base_url=api_base_url
        )
        print("✓ Success")
        
        # Database information
        info = agent.get_feature_info()
        print("\n📊 KNOWLEDGE DATABASE")
        print("─" * 50)
        print(f"• Total Features Available: {info['total_features']}")
        print(f"• Sample Dataset Size: {info['sample_count']} patients")
        print(f"• Feature Categories: {len(info['categories'])} types")
        print()
        
        for category, count in info['categories'].items():
            print(f"  ├─ {category:<20}: {count:3d} features")
        print()
        
        # Professional demo
        print("🔬 DEMONSTRATION ANALYSIS")
        print("═" * 70)
        
        # Question
        example_question = "What does Ripley_K_mean_Habitat08 mean?"
        print("\n📋 RESEARCH QUESTION:")
        print("┌" + "─" * 68 + "┐")
        print(f"│ {example_question:<66} │")
        print("└" + "─" * 68 + "┘")
        
        # Processing
        print("\n🤖 EXPERT ANALYSIS:")
        print("─" * 70)
        print("⏳ Processing query through CANVAS knowledge base...")
        
        # Answer
        result = agent.ask_question(example_question)
        print("\n📄 EXPERT RESPONSE:")
        print("─" * 70)
        print(result)
        print("─" * 70)
        
        # Interactive mode invitation
        print("\n🎯 INTERACTIVE MODE")
        print("═" * 70)
        print("Ready to answer your research questions!")
        print("Type any question about spatial features, habitats, or immunotherapy analysis.")
        print("Commands: 'help' for assistance, 'quit' to exit")
        print()
        
        agent.interactive_mode()
        
    except Exception as e:
        print("✗ Failed")
        print(f"\n❌ ERROR: {e}")
        print("Please check your configuration and try again.")

if __name__ == "__main__":
    quick_demo()
