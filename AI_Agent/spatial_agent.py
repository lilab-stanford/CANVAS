import json
from openai import OpenAI
from concurrent.futures import ThreadPoolExecutor, wait
from functools import partial
from tqdm import tqdm
import os
import argparse
import pandas as pd
import time
from typing import Dict, List, Any, Tuple
import docx

class SpatialFeatureAgent:
    """
    A GPT-based agent for interpreting CANVAS-derived spatial features
    """
    
    def __init__(self, api_key=None, model_name="gpt-4o", temperature=0.3, 
                 num_workers=8, timeout_duration=60, retry_attempts=3, api_base_url=None):
        """
        Initialize the Spatial Feature Agent
        
        Args:
            api_key: OpenAI API key (required if OPENAI_API_KEY env var not set)
            model_name: GPT model to use
            temperature: Temperature for generation
            num_workers: Number of parallel workers
            timeout_duration: Request timeout in seconds
            retry_attempts: Number of retry attempts
            api_base_url: Optional custom API base URL; omit when using the official OpenAI endpoint
        """
        self.api_key = api_key or os.getenv('OPENAI_API_KEY')
        if not self.api_key:
            raise ValueError(
                "OpenAI API key is required. Provide it via the api_key argument or set the OPENAI_API_KEY environment variable."
            )

        # Allow optional override of API base URL for non-official deployments
        self.api_base_url = api_base_url or os.getenv('OPENAI_BASE_URL') or os.getenv('OPENAI_API_BASE')

        client_kwargs = {"api_key": self.api_key}
        if self.api_base_url:
            client_kwargs["base_url"] = self.api_base_url

        self.client = OpenAI(**client_kwargs)
            
        self.model_name = model_name
        self.temperature = temperature
        self.num_workers = num_workers
        self.timeout_duration = timeout_duration
        self.retry_attempts = retry_attempts
        self.miss_index = []
        
        # Load context data and create system prompt
        self.context_data = self._load_context_data()
        self.system_prompt = self._create_system_prompt()
        
    def _load_context_data(self) -> Dict[str, Any]:
        """Load all context data from the data directory"""
        context = {}
        
        try:
            # Load feature information from Feature_annotation.xlsx
            # Sheet1: Feature type definitions
            feature_types_df = pd.read_excel('data/Feature_annotation.xlsx', sheet_name='Sheet1')
            # Sheet2: All features with their types
            features_df = pd.read_excel('data/Feature_annotation.xlsx', sheet_name='Sheet2')
            
            # Extract feature list from Sheet2
            context['available_features'] = features_df['Feature_ID'].tolist()
            context['total_features'] = len(features_df)
            context['sample_count'] = 149  # Fixed sample count
            
            # Create feature categories based on Feature_Type from Sheet2
            feature_type_mapping = dict(zip(features_df['Feature_ID'], features_df['Feature_Type']))
            context['feature_categories'] = {}
            
            # Group features by their types
            for feature_id, feature_type in feature_type_mapping.items():
                if feature_type not in context['feature_categories']:
                    context['feature_categories'][feature_type] = []
                context['feature_categories'][feature_type].append(feature_id)
            
            # Load feature type definitions from Sheet1
            try:
                feature_definitions = {}
                for i, row in feature_types_df.iterrows():
                    if i == 0:  # Skip header row
                        continue
                    if i <= 6:  # Process the 6 feature types
                        feature_type = row.iloc[1] if pd.notna(row.iloc[1]) else ""
                        definition = row.iloc[2] if pd.notna(row.iloc[2]) else ""
                        biological = row.iloc[3] if pd.notna(row.iloc[3]) else ""
                        
                        if feature_type:
                            feature_definitions[feature_type] = {
                                'definition': definition,
                                'biological_interpretation': biological
                            }
                
                context['feature_type_definitions'] = feature_definitions
            except Exception as e:
                print(f"Warning: Could not load feature type definitions: {e}")
                context['feature_type_definitions'] = {}
            
            # Load habitat annotations from docx
            try:
                doc = docx.Document('data/Habitat_annotation.docx')
                habitat_text = '\n'.join([p.text for p in doc.paragraphs if p.text.strip()])
                context['habitat_annotations'] = habitat_text
            except Exception as e:
                print(f"Warning: Could not load habitat annotations: {e}")
                context['habitat_annotations'] = ""
                
        except Exception as e:
            print(f"Warning: Could not load feature annotation data: {e}")
            context['available_features'] = []
            context['feature_categories'] = {}
        
        return context
    
    def _create_system_prompt(self) -> str:
        """Create the system prompt based on the instruction document"""
        system_prompt = """
🔬 You are an expert AI consultant specializing in CANVAS-derived spatial features for NSCLC immunotherapy analysis.

🎯 Your Role:
You can answer ANY question related to:
- Spatial feature interpretation and analysis
- Habitat biology and cellular composition
- Immunotherapy response prediction
- Spatial dispersion and statistical methods
- Clinical implications and prognostic significance
- Comparative analysis between features/habitats
- General questions about tumor microenvironment

🧠 When interpreting specific features, provide structured analysis with these 5 dimensions:

1. Category
Identify which of the six major domains the feature belongs to:
- Composition
- Diversity  
- Spatial dispersion
- Interaction
- Distance
- Transition

2. Cellular Composition of Associated Habitat
Describe the dominant or enriched cell types in the corresponding habitat, referencing known immune or stromal populations:
- Tcyto, Th, B cell, Plasma cell → immune-activating
- CAF, Neutrophil, M2 → immunosuppressive  
- EC, DCs, Mono → context-dependent

3. Spatial Property Description
Explain what the spatial feature measures and what a high or low value indicates. For example:
- Ripley_K_mean: degree of spatial clustering or dispersion across radii
- Kernel density: local concentration
- J_function: spatial regularity
- Pairwise interaction: proximity-based coupling across habitats
- STE: entropy-based intermixing

4. Topological Coupling Tendency
Comment on which other habitats this feature tends to correlate or co-occur with, based on Jaccard proximity, interaction networks, or shared prognostic modules (e.g., H03–H08 co-localization, H01–H04–H06 exclusion).

5. Biological and Clinical Implication
Synthesize its potential role in:
- Immune activation vs. suppression
- Prognostic outcome (protective or risk)
- Predictive relevance to immunotherapy response

🏠 Habitat Biological Annotations:
- H01 (Habitat01): Tumorogenic Core - malignant cells, tumor-intrinsic features
- H02 (Habitat02): Macrophage Enriched - M1/M2 macrophages, immune regulation
- H03 (Habitat03): B-cell Enriched - B cells, humoral immunity, potential TLS formation
- H04 (Habitat04): Fibrotic Activity Hub - CAFs, ECM remodeling, immunosuppressive stroma
- H05 (Habitat05): Plasma cell Enriched - plasma cells, antibody production
- H06 (Habitat06): Neutrophil Prominent - neutrophils, often immunosuppressive N2-like TANs
- H07 (Habitat07): Tumor Interface - tumor-stroma boundary, transitional zone
- H08 (Habitat08): T-lymphonic Enriched - CD4+/CD8+ T cells, cellular immunity, TLS
- H09 (Habitat09): Pan-immune Active Zone - diverse immune cells, immune activation hub
- H10 (Habitat10): Vasculature Niche - endothelial cells, angiogenesis, vascular features

📊 Feature Categories Knowledge:

1. **Composition (10 features)**: frequency_Habitat01-10
   - Definition: Quantifies the frequency or proportion of each habitat per image
   - Biology: Reflects abundance of specific microenvironmental states; e.g., high H01 implies tumor core dominance
   - Calculation: Proportion of habitat label across all patches per sample

2. **Diversity (6 features)**: div_Richness, div_Shannon, div_Simpson, div_Inv_Simpson, div_Pielou, div_Fisher_alpha
   - Definition: Measures complexity of habitat composition (e.g., Shannon, Simpson, Richness, etc.)
   - Biology: High diversity (e.g., high Shannon with enrichment of H02/H03/H08) is linked to better immune surveillance and broader immune repertoire
   - Calculation: vegan::diversity, spehabitatumber, fisher.alpha, Shannon = entropy of habitats

3. **Spatial Dispersion (~90 features)**: Ripley_K/L_mean, G_mean, F_mean, J_mean, Pair_corr_g_mean, Clark_Evans, Quadrat_chisq, Kernel_density_mean
   - Definition: Distribution within habitats (Ripley's K/L, Clark-Evans, F/G/J, Kernel, Quadrat)
   - Biology: Aggregation (high K/L) = immune barrier; uniformity (Clark-Evans) = favorable environment
   - Calculation: Ripley's K/L, Clark-Evans, Quadrat, Kernel density via point pattern analysis

4. **Interaction (100 features)**: cci_HabitatX_HabitatY
   - Definition: Captures cell-cell neighborhood relationships (e.g., GNN edge weight, SCIMAP Z-score)
   - Biology: High interaction scores (e.g., H01–H03 or H01–H08) may indicate immune cell clustering or immune surveillance–associated niches
   - Calculation: SCIMAP spatial_interaction, Z-score calculation, or GNN-derived edge weights

5. **Distance (55 features)**: dis_HabitatX_HabitatY
   - Definition: Pairwise habitat distances based on nearest neighbor or exclusion patterns
   - Biology: Shorter distances suggest active immune infiltration, whereas longer distances may reflect immune exclusion or spatial segregation
   - Calculation: Patch-level pairwise habitat average distances, then ranking

6. **Transition (1 feature)**: SpatialTransitionEntropy
   - Definition: Spatial Transition Entropy between habitats across image patches
   - Biology: High entropy indicates spatial complexity and mixed microenvironments; low entropy suggests structured or compartmentalized tissue
   - Calculation: Entropy of patch-wise habitat transition matrix (non-distance-based entropy)

🔬 Known Clinical Associations:
- Protective features (β < 0): H03, H05, H07, H08 proximity; organized immune zones; TLS formation
- Risk features (β > 0): H01, H04, H06 clustering; immunosuppressive configurations; immune exclusion
- H03-H08 co-localization → coordinated adaptive immune response → favorable outcomes
- H01-H04-H06 clustering → immunosuppressive TME → poor prognosis
- High spatial diversity → organized immune compartmentalization → ICB responsiveness

💡 You can answer questions like:
- "What does Ripley_K_mean_Habitat08 mean?" (structured 5-dimension analysis)
- "Compare H03 and H08 habitats" (comparative analysis)
- "Which features predict good immunotherapy response?" (clinical insights)
- "Explain spatial dispersion in general" (educational content)
- "What's the difference between Ripley K and L functions?" (methodology)
- "Show me interaction features for Habitat01" (data queries)

Always provide concrete, biologically grounded responses. Be conversational and educational while maintaining scientific accuracy.
"""
        return system_prompt.strip()

    def _process_single_query(self, message_data: Tuple[int, str]) -> Tuple[int, str]:
        """Process single query request"""
        index, user_query = message_data
        try:
            completion = self.client.chat.completions.create(
                model=self.model_name,
                messages=[
                    {"role": "system", "content": self.system_prompt},
                    {"role": "user", "content": user_query}
                ],
                temperature=self.temperature,
                max_tokens=1500
            )
            return (index, completion.choices[0].message.content.strip())
        except Exception as e:
            print(f"Error processing query '{user_query}': {e}")
            self.miss_index.append(index)
            return (index, None)
    
    def _process_query_batch(self, query_list: List[str]) -> List[Tuple[int, str]]:
        """Process multiple query requests"""
        results = []
        executor = ThreadPoolExecutor(max_workers=self.num_workers)
        
        # Create indexed list for processing
        indexed_queries = [(index, query) for index, query in enumerate(query_list)]
        
        try:
            # Submit all tasks
            future_to_query = {
                executor.submit(self._process_single_query, query_data): query_data 
                for query_data in indexed_queries
            }
            
            # Process results with timeout and retry
            for _ in range(self.retry_attempts):
                done, not_done = wait(future_to_query.keys(), timeout=self.timeout_duration)
                
                # Cancel timed out futures
                for future in not_done:
                    future.cancel()
                
                # Collect results
                for future in done:
                    if future.done() and not future.cancelled():
                        try:
                            result = future.result()
                            results.append(result)
                        except Exception as e:
                            print(f"Error getting future result: {e}")
                
                if len(not_done) == 0:
                    break
                    
                # Retry failed requests
                print(f"Retrying {len(not_done)} failed requests...")
                future_to_query = {
                    executor.submit(self._process_single_query, future_to_query[future]): 
                    future_to_query[future] for future in not_done
                }
                
        except Exception as e:
            print(f"Error in batch processing: {e}")
        finally:
            executor.shutdown(wait=True)
            
        return results
    
    def ask_question(self, question: str) -> str:
        """
        Ask any question about spatial features, habitats, or immunotherapy analysis
        
        Args:
            question: Any question related to spatial analysis
            
        Returns:
            Expert response to the question
        """
        results = self._process_query_batch([question])
        if results and results[0][1] is not None:
            return results[0][1]
        else:
            return f"Failed to process question: {question}"
    
    def interpret_feature(self, feature_name: str) -> str:
        """
        Interpret a single spatial feature (legacy method for backward compatibility)
        
        Args:
            feature_name: Name of the spatial feature to interpret
            
        Returns:
            Structured interpretation of the feature
        """
        if feature_name not in self.context_data.get('available_features', []):
            return f"Warning: '{feature_name}' is not found in the 262 available features. Please check the feature name."
        
        question = f"Please provide a detailed interpretation of the spatial feature: {feature_name}"
        return self.ask_question(question)
    
    def interpret_features_batch(self, feature_names: List[str], output_file: str = None) -> Dict[str, str]:
        """
        Interpret multiple spatial features in batch
        
        Args:
            feature_names: List of feature names to interpret
            output_file: Optional output file to save results
            
        Returns:
            Dictionary mapping feature names to interpretations
        """
        # Validate features
        valid_features = []
        invalid_features = []
        
        for feature in feature_names:
            if feature in self.context_data.get('available_features', []):
                valid_features.append(feature)
            else:
                invalid_features.append(feature)
        
        if invalid_features:
            print(f"Warning: The following features are not valid: {invalid_features}")
        
        if not valid_features:
            return {}
        
        print(f"Processing {len(valid_features)} features...")
        
        # Create questions for batch processing
        questions = [f"Please provide a detailed interpretation of the spatial feature: {feature}" for feature in valid_features]
        
        # Process questions in batch
        results = self._process_query_batch(questions)
        
        # Sort results by index
        results.sort(key=lambda x: x[0])
        
        # Create output dictionary
        interpretation_dict = {}
        for i, (index, interpretation) in enumerate(results):
            if i < len(valid_features):
                feature_name = valid_features[index]
                interpretation_dict[feature_name] = interpretation if interpretation else f"Failed to interpret: {feature_name}"
        
        # Handle missing results
        for i, feature in enumerate(valid_features):
            if feature not in interpretation_dict:
                interpretation_dict[feature] = f"Failed to interpret: {feature}"
        
        # Save to file if specified
        if output_file:
            self._save_results(interpretation_dict, output_file)
        
        print(f"Successfully interpreted {len([v for v in interpretation_dict.values() if not v.startswith('Failed')])} out of {len(valid_features)} features")
        
        return interpretation_dict
    
    def _save_results(self, results: Dict[str, str], output_file: str):
        """Save interpretation results to file"""
        try:
            # Create output directory if it doesn't exist
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            
            if output_file.endswith('.json'):
                with open(output_file, 'w', encoding='utf-8') as f:
                    json.dump(results, f, indent=2, ensure_ascii=False)
            elif output_file.endswith('.csv'):
                df = pd.DataFrame([
                    {'Feature': feature, 'Interpretation': interpretation}
                    for feature, interpretation in results.items()
                ])
                df.to_csv(output_file, index=False)
            else:
                # Save as text file
                with open(output_file, 'w', encoding='utf-8') as f:
                    for feature, interpretation in results.items():
                        f.write(f"=== {feature} ===\n")
                        f.write(f"{interpretation}\n\n")
            
            print(f"Results saved to: {output_file}")
        except Exception as e:
            print(f"Error saving results: {e}")
    
    def list_available_features(self, category: str = None) -> List[str]:
        """List available features, optionally filtered by category"""
        if not category:
            return self.context_data.get('available_features', [])
        
        return self.context_data.get('feature_categories', {}).get(category, [])
    
    def get_feature_info(self) -> Dict[str, Any]:
        """Get summary information about available features"""
        return {
            'total_features': self.context_data.get('total_features', 0),
            'sample_count': self.context_data.get('sample_count', 0),
            'categories': {k: len(v) for k, v in self.context_data.get('feature_categories', {}).items()}
        }
    
    def interactive_mode(self):
        """Start interactive mode for open-ended questions"""
        # Show some example questions
        examples = [
            "What does Ripley_K_mean_Habitat08 mean?",
            "Compare H03 and H08 habitats", 
            "Which features predict immunotherapy response?",
            "Explain Ripley K function",
            "What are interaction features?",
            "Show me all diversity features"
        ]
        
        print("💡 EXAMPLE RESEARCH QUESTIONS:")
        print("─" * 50)
        for i, example in enumerate(examples, 1):
            print(f"  {i}. {example}")
        print("\n" + "═" * 70)
        
        session_count = 0
        while True:
            try:
                # Professional question prompt
                print(f"\n🔬 QUERY #{session_count + 1}")
                print("─" * 30)
                user_input = input("❓ Research Question: ").strip()
                
                if user_input.lower() in ['quit', 'exit', 'q']:
                    print("\n🔚 SESSION TERMINATED")
                    print("─" * 40)
                    print("Thank you for using CANVAS Expert!")
                    print("─" * 40)
                    break
                elif user_input.lower() == 'help':
                    self._show_help()
                elif user_input.lower() == 'info':
                    info = self.get_feature_info()
                    print("\n📊 DATABASE SUMMARY:")
                    print("─" * 40)
                    print(f"• Total Features: {info['total_features']}")
                    print(f"• Patient Samples: {info['sample_count']}")
                    print(f"• Feature Categories: {len(info['categories'])}")
                    print("\n📋 CATEGORY BREAKDOWN:")
                    for category, count in info['categories'].items():
                        print(f"  ├─ {category:18}: {count:3d} features")
                elif user_input.lower().startswith('list'):
                    parts = user_input.split()
                    category = parts[1] if len(parts) > 1 else None
                    features = self.list_available_features(category)
                    
                    print(f"\n📋 FEATURE CATALOG ({len(features)} items):")
                    print("─" * 50)
                    for i, feature in enumerate(features[:20]):  # Show first 20
                        print(f"  {i+1:2d}. {feature}")
                    if len(features) > 20:
                        print(f"\n     ... plus {len(features)-20} additional features")
                        print("     Use 'list [category]' to filter by category")
                elif user_input:
                    session_count += 1
                    
                    # Processing indication
                    print(f"\n🤖 ANALYSIS IN PROGRESS (Query #{session_count})")
                    print("─" * 60)
                    print("⏳ Processing through CANVAS knowledge base...")
                    
                    # Get response
                    response = self.ask_question(user_input)
                    
                    # Clean, professional response format
                    print(f"\n📋 RESEARCH QUESTION:")
                    print(f"   {user_input}")
                    print(f"\n📄 EXPERT RESPONSE:")
                    print("─" * 60)
                    print(response)
                    print("─" * 60)
                    print("✅ Analysis complete. Ready for next query.\n")
                
            except KeyboardInterrupt:
                print("\n\n⏹️  SESSION INTERRUPTED")
                print("─" * 40)
                print("Session ended by user")
                break
            except Exception as e:
                print(f"\n❌ SYSTEM ERROR: {e}")
                print("Please try rephrasing your question or contact support.")
    
    def _show_help(self):
        """Show help information"""
        print("\n📚 HELP DOCUMENTATION")
        print("═" * 50)
        
        print("\n📋 SYSTEM COMMANDS:")
        print("─" * 40)
        commands = [
            ("help", "Display this help documentation"),
            ("info", "Show database statistics and overview"),
            ("list", "Display all features (first 20 items)"),
            ("list [category]", "Filter features by category"),
            ("quit/exit/q", "Terminate the analysis session")
        ]
        
        for cmd, desc in commands:
            print(f"  ├─ {cmd:<15} : {desc}")
        
        print("\n🔬 RESEARCH QUESTION CATEGORIES:")
        print("─" * 50)
        
        categories = [
            ("Feature Analysis", [
                "What does Ripley_K_mean_Habitat08 mean?",
                "Interpret frequency_Habitat03",
                "Explain div_Shannon diversity index"
            ]),
            ("Comparative Analysis", [
                "Compare H03 and H08 habitats",
                "Difference between Ripley K and L functions?",
                "Which is better: Shannon or Simpson diversity?"
            ]),
            ("Clinical Insights", [
                "Which features predict immunotherapy response?",
                "What makes a tumor immunosuppressive?",
                "How do T cells and B cells interact spatially?"
            ]),
            ("Data Queries", [
                "Show me all interaction features for Habitat01",
                "List diversity features",
                "What spatial dispersion are available?"
            ]),
            ("Educational Content", [
                "Explain spatial point pattern analysis",
                "What are tertiary lymphoid structures?",
                "How does tumor microenvironment affect treatment?"
            ])
        ]
        
        for category, examples in categories:
            print(f"\n🎯 {category}:")
            for example in examples:
                print(f"    • {example}")
        
        print("\n" + "═" * 50)
        print("💡 TIP: Ask questions in natural language - no special syntax required!")
        print("═" * 50)

def parse_args():
    parser = argparse.ArgumentParser(description='Spatial Feature Interpretation Agent')
    parser.add_argument('--api_key', type=str, help='OpenAI API key')
    parser.add_argument('--api_url', type=str, help='Custom API base URL')
    parser.add_argument('--model', type=str, default="gpt-4o", help='GPT model to use')
    parser.add_argument('--feature', type=str, help='Single feature to interpret')
    parser.add_argument('--features_file', type=str, help='File containing list of features to interpret')
    parser.add_argument('--output', type=str, help='Output file for batch results')
    parser.add_argument('--interactive', action='store_true', default=True, help='Start interactive mode')
    parser.add_argument('--num_workers', type=int, default=8, help='Number of parallel workers')
    parser.add_argument('--timeout', type=int, default=60, help='Request timeout in seconds')
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Initialize agent
    try:
        agent = SpatialFeatureAgent(
            api_key=args.api_key,
            model_name=args.model,
            api_base_url=args.api_url,
            num_workers=args.num_workers,
            timeout_duration=args.timeout
        )
    except ValueError as exc:
        print(f"Configuration error: {exc}")
        return
    
    if args.interactive:
        # Start interactive mode
        agent.interactive_mode()
    elif args.feature:
        # Interpret single feature
        print(f"Interpreting feature: {args.feature}")
        result = agent.interpret_feature(args.feature)
        print(f"\nInterpretation:\n{result}")
    elif args.features_file:
        # Batch interpret features from file
        try:
            with open(args.features_file, 'r') as f:
                features = [line.strip() for line in f if line.strip()]
            
            print(f"Loading {len(features)} features from {args.features_file}")
            results = agent.interpret_features_batch(features, args.output)
            
            if not args.output:
                for feature, interpretation in results.items():
                    print(f"\n=== {feature} ===")
                    print(interpretation)
                    print("-" * 80)
                    
        except Exception as e:
            print(f"Error processing features file: {e}")
    else:
        # Show feature info and start interactive mode
        info = agent.get_feature_info()
        print("=== Spatial Feature Agent ===")
        print(f"Total features available: {info['total_features']}")
        print("Feature categories:")
        for category, count in info['categories'].items():
            print(f"  {category}: {count}")
        print("\nUse --interactive to start interactive mode")
        print("Use --feature <feature_name> to interpret a single feature")
        print("Use --help for more options")

if __name__ == "__main__":
    main() 
