# Build a content builder agent

[Tutorials](/oss/python/deepagents/data-analysis)[Deep Agents](/oss/python/deepagents/data-analysis)

# Build a content builder agent
Copy page

Build a content writing agent with brand memory, skills, subagents, and image generationCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

This guide demonstrates how to build a content writing agent from scratch using [Deep Agents](/oss/python/deepagents).
The agent you build will:

{' ' * (self.list_depth - 1)}- Load voice and workflow rules from `AGENTS.md` and skill folders

{' ' * (self.list_depth - 1)}- Delegate web research to a specialized subagent with `web_search`

{' ' * (self.list_depth - 1)}- Draft blog or social content following the loaded skill

{' ' * (self.list_depth - 1)}- Generate cover or social images with Gemini and save files under the project directory

The code in this tutorial wires in image generation tools and a filesystem backend so the agent can read and write posts, research notes, and images under the project directory. For the full runnable project, see the [content-builder-agent](https://github.com/langchain-ai/deepagents/tree/main/examples/content-builder-agent) example.

### [​](#key-concepts)Key concepts

This tutorial covers:

{' ' * (self.list_depth - 1)}- [Long-term memory](/oss/python/deepagents/long-term-memory) for TODO

{' ' * (self.list_depth - 1)}- [Skills](/oss/python/deepagents/skills) for TODO

{' ' * (self.list_depth - 1)}- [Subagents](/oss/python/deepagents/subagents) for TODO

{' ' * (self.list_depth - 1)}- [Filesystem backends](/oss/python/deepagents/backends) for file read and write

{' ' * (self.list_depth - 1)}- Custom [tools](/oss/python/langchain/tools) for search and image generation

## [​](#prerequisites)Prerequisites

API keys:

{' ' * (self.list_depth - 1)}- Anthropic (Claude) or other provider API key

{' ' * (self.list_depth - 1)}- Google (Gemini) for image generation with `gemini-2.5-flash-image`

{' ' * (self.list_depth - 1)}- [Tavily](https://www.tavily.com/) for web search (free tier)

{' ' * (self.list_depth - 1)}- [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-deepagents-content-builder) for tracing (optional)

Python 3.11 or later.

## [​](#setup)Setup

1[](#)

Create project directory
```
mkdir content-builder-agent
cd content-builder-agent

```
2[](#)

Install dependenciespipuv
```
pip install deepagents google-genai pillow pyyaml rich tavily-python langchain

```
Pin `deepagents` to a supported range in your own project (for example `>=0.3.5,<0.4.0`) to match the upstream example.3[](#)

Set API keys
```
export ANTHROPIC_API_KEY="your_anthropic_api_key"
export GOOGLE_API_KEY="your_google_api_key"
export TAVILY_API_KEY="your_tavily_api_key" # Optional
export LANGSMITH_API_KEY="your_langsmith_api_key" # Optional

```

## [​](#add-configuration-files)Add configuration files

The example keeps behavior in three kinds of files: memory, skills, and subagent definitions.
1[](#)

Add AGENTS.mdCreate `AGENTS.md` in the project root.
When you later create the agent and specify this file as part of the [memory](/oss/python/deepagents/long-term-memory) parameter, it gets loaded this into the system prompt so brand voice and research expectations apply to every run.
```
# Content Writer Agent

You are a content writer for a technology company. Your job is to create engaging, informative content that educates readers about AI, software development, and emerging technologies.

## Brand Voice

- **Professional but approachable**: Write like a knowledgeable colleague, not a textbook
- **Clear and direct**: Avoid jargon unless necessary; explain technical concepts simply
- **Confident but not arrogant**: Share expertise without being condescending
- **Engaging**: Use concrete examples, analogies, and stories to illustrate points

## Writing Standards

1. **Use active voice**: "The agent processes requests" not "Requests are processed by the agent"
2. **Lead with value**: Start with what matters to the reader, not background
3. **One idea per paragraph**: Keep paragraphs focused and scannable
4. **Concrete over abstract**: Use specific examples, numbers, and case studies
5. **End with action**: Every piece should leave the reader knowing what to do next

## Content Pillars

Our content focuses on:
- AI agents and automation
- Developer tools and productivity
- Software architecture and best practices
- Emerging technologies and trends

## Formatting Guidelines

- Use headers (H2, H3) to break up long content
- Include code examples where relevant (with syntax highlighting)
- Add bullet points for lists of 3+ items
- Keep sentences under 25 words when possible
- Include a clear call-to-action at the end

## Research Requirements

Before writing on any topic:
1. Use the `researcher` subagent for in-depth topic research
2. Gather at least 3 credible sources
3. Identify the key points readers need to understand
4. Find concrete examples or case studies to illustrate concepts

```
See all 42 linesTo make this agent comply with your own tone, pillars, and formatting rules, update the text in `AGENTS.md`.2[](#)

Add subagents.yamlCreate a file called `subagents.yaml`.
Then add the following txt which describes a `researcher` subagent with a Tavily-backed `web_search` tool, a Haiku model id, and instructions to save findings to paths you specify when delegating from the main agent:
```
# Subagent definitions
# These are loaded by content_writer.py and wired up with tools

researcher:
 description: >
 ALWAYS use this first to research any topic before writing content.
 Searches the web for current information, statistics, and sources.
 When delegating, tell it the topic AND the file path to save results
 (e.g., 'Research renewable energy and save to research/renewable-energy.md').
 model: anthropic:claude-haiku-4-5-20251001
 system_prompt: |
 You are a research assistant. You have access to web_search and write_file tools.

 ## Your Tools
 - web_search(query, max_results=5, topic="general") - Search the web
 - write_file(file_path, content) - Save your findings

 ## Your Process
 1. Use web_search to find information on the topic
 2. Make 2-3 targeted searches with specific queries
 3. Gather key statistics, quotes, and examples
 4. Save findings to the file path specified in your task

 ## Important
 - The user will tell you WHERE to save the file - use that exact path
 - Always include source URLs in your findings
 - Keep findings concise but informative
 tools:
 - web_search

```
See all 29 linesThe file gets passed as an argument later when creating the deep agent.3[](#)

Add skillsCreate a `skills/` directory. Each skill is a folder containing a `SKILL.md` file with YAML frontmatter (`name`, `description`) and instructions for the skill.Create `skills/blog-post/SKILL.md` and copy the following text into it which contains information on crating long-form posts, optimizing content for SEO, and generating cover images.
```
---
name: blog-post
description: Writes and structures long-form blog posts, creates tutorial outlines, and optimizes content for SEO with cover image generation. Use when the user asks to write a blog post, article, how-to guide, tutorial, technical writeup, thought leadership piece, or long-form content.
---

# Blog Post Writing Skill

## Research First (Required)

**Before writing any blog post, you MUST delegate research:**

1. Use the `task` tool with `subagent_type: "researcher"`
2. In the description, specify BOTH the topic AND where to save:

```
task(
 subagent_type="researcher",
 description="Research [TOPIC]. Save findings to research/[slug].md"
)
```

Example:
```
task(
 subagent_type="researcher",
 description="Research the current state of AI agents in 2025. Save findings to research/ai-agents-2025.md"
)
```

3. After research completes, read the findings file before writing

## Output Structure (Required)

**Every blog post MUST have both a post AND a cover image:**

```
blogs/
└── <slug>/
 ├── post.md # The blog post content
 └── hero.png # REQUIRED: Generated cover image
```

Example: A post about "AI Agents in 2025" → `blogs/ai-agents-2025/`

**You MUST complete both steps:**
1. Write the post to `blogs/<slug>/post.md`
2. Generate a cover image using `generate_image` and save to `blogs/<slug>/hero.png`

**A blog post is NOT complete without its cover image.**

## Blog Post Structure

Every blog post should follow this structure:

### 1. Hook (Opening)
- Start with a compelling question, statistic, or statement
- Make the reader want to continue
- Keep it to 2-3 sentences

### 2. Context (The Problem)
- Explain why this topic matters
- Describe the problem or opportunity
- Connect to the reader's experience

### 3. Main Content (The Solution)
- Break into 3-5 main sections with H2 headers
- Each section covers one key point
- Include code examples, diagrams, or screenshots where helpful
- Use bullet points for lists

### 4. Practical Application
- Show how to apply the concepts
- Include step-by-step instructions if applicable
- Provide code snippets or templates

### 5. Conclusion & CTA
- Summarize key takeaways (3 bullets max)
- End with a clear call-to-action
- Link to related resources

## Cover Image Generation

After writing the post, generate a cover image using the `generate_cover` tool:

```
generate_cover(prompt="A detailed description of the image...", slug="your-blog-slug")
```

The tool saves the image to `blogs/<slug>/hero.png`.

### Writing Effective Image Prompts

Structure your prompt with these elements:

1. **Subject**: What is the main focus? Be specific and concrete.
2. **Style**: Art direction (minimalist, isometric, flat design, 3D render, watercolor, etc.)
3. **Composition**: How elements are arranged (centered, rule of thirds, symmetrical)
4. **Color palette**: Specific colors or mood (warm earth tones, cool blues and purples, high contrast)
5. **Lighting/Atmosphere**: Soft diffused light, dramatic shadows, golden hour, neon glow
6. **Technical details**: Aspect ratio considerations, negative space for text overlay

### Example Prompts

**For a technical blog post:**
```
Isometric 3D illustration of interconnected glowing cubes representing AI agents, each cube has subtle circuit patterns. Cubes connected by luminous data streams. Deep navy background (#0a192f) with electric blue (#64ffda) and soft purple (#c792ea) accents. Clean minimal style, lots of negative space at top for title. Professional tech aesthetic.
```

**For a tutorial/how-to:**
```
Clean flat illustration of hands typing on a keyboard with abstract code symbols floating upward, transforming into lightbulbs and gears. Warm gradient background from soft coral to light peach. Friendly, approachable style. Centered composition with space for text overlay.
```

**For thought leadership:**
```
Abstract visualization of a human silhouette profile merging with geometric neural network patterns. Split composition - organic watercolor texture on left transitioning to clean vector lines on right. Muted sage green and warm terracotta color scheme. Contemplative, forward-thinking mood.
```

## SEO Considerations

- Include the main keyword in the title and first paragraph
- Use the keyword naturally 3-5 times throughout
- Keep the title under 60 characters
- Write a meta description (150-160 characters)

## Quality Checklist

Before finishing:
- [ ] Post saved to `blogs/<slug>/post.md`
- [ ] Hero image generated at `blogs/<slug>/hero.png`
- [ ] Hook grabs attention in first 2 sentences
- [ ] Each section has a clear purpose
- [ ] Conclusion summarizes key points
- [ ] CTA tells reader what to do next

```
See all 134 linesNext, create `skills/social-media/SKILL.md` and copy the following text into it which contains information on drafting social media posts and generating accompanying imagery:
```
---
name: social-media
description: Drafts engaging social media posts, writes hooks, suggests hashtags, creates thread structures, and generates companion images. Use when the user asks to write a LinkedIn post, tweet, Twitter/X thread, social media caption, social post, or repurpose content for social platforms.
---

# Social Media Content Skill

## Research First (Required)

**Before writing any social media content, you MUST delegate research:**

1. Use the `task` tool with `subagent_type: "researcher"`
2. In the description, specify BOTH the topic AND where to save:

```
task(
 subagent_type="researcher",
 description="Research [TOPIC]. Save findings to research/[slug].md"
)
```

Example:
```
task(
 subagent_type="researcher",
 description="Research renewable energy trends in 2025. Save findings to research/renewable-energy.md"
)
```

3. After research completes, read the findings file before writing

## Output Structure (Required)

**Every social media post MUST have both content AND an image:**

**LinkedIn posts:**
```
linkedin/
└── <slug>/
 ├── post.md # The post content
 └── image.png # REQUIRED: Generated visual
```

**Twitter/X threads:**
```
tweets/
└── <slug>/
 ├── thread.md # The thread content
 └── image.png # REQUIRED: Generated visual
```

Example: A LinkedIn post about "prompt engineering" → `linkedin/prompt-engineering/`

**You MUST complete both steps:**
1. Write the content to the appropriate path
2. Generate an image using `generate_image` and save alongside the post

**A social media post is NOT complete without its image.**

## Platform Guidelines

### LinkedIn

**Format:**
- 1,300 character limit (show more after ~210 chars)
- First line is crucial - make it hook
- Use line breaks for readability
- 3-5 hashtags at the end

**Tone:**
- Professional but personal
- Share insights and learnings
- Ask questions to drive engagement
- Use "I" and share experiences

**Structure:**
```
[Hook - 1 compelling line]

[Empty line]

[Context - why this matters]

[Empty line]

[Main insight - 2-3 short paragraphs]

[Empty line]

[Call to action or question]

#hashtag1 #hashtag2 #hashtag3
```

### Twitter/X

**Format:**
- 280 character limit per tweet
- Threads for longer content (use 1/🧵 format)
- No more than 2 hashtags per tweet

**Thread Structure:**
```
1/🧵 [Hook - the main insight]

2/ [Supporting point 1]

3/ [Supporting point 2]

4/ [Example or evidence]

5/ [Conclusion + CTA]
```

## Image Generation

Every social media post needs an eye-catching image. Use the `generate_social_image` tool:

```
generate_social_image(prompt="A detailed description...", platform="linkedin", slug="your-post-slug")
```

The tool saves the image to `<platform>/<slug>/image.png`.

### Social Image Best Practices

Social images need to work at small sizes in crowded feeds:
- **Bold, simple compositions** - one clear focal point
- **High contrast** - stands out when scrolling
- **No text in image** - too small to read, platforms add their own
- **Square or 4:5 ratio** - works across platforms

### Writing Effective Prompts

Include these elements:

1. **Single focal point**: One clear subject, not a busy scene
2. **Bold style**: Vibrant colors, strong shapes, high contrast
3. **Simple background**: Solid color, gradient, or subtle texture
4. **Mood/energy**: Match the post tone (inspiring, urgent, thoughtful)

### Example Prompts

**For an insight/tip post:**
```
Single glowing lightbulb floating against a deep purple gradient background, lightbulb made of interconnected golden geometric lines, rays of soft light emanating outward. Minimal, striking, high contrast. Square composition.
```

**For announcements/news:**
```
Abstract rocket ship made of colorful geometric shapes launching upward with a trail of particles. Bright coral and teal color scheme against clean white background. Energetic, celebratory mood. Bold flat illustration style.
```

**For thought-provoking content:**
```
Two overlapping translucent circles, one blue one orange, creating a glowing intersection in the center. Represents collaboration or intersection of ideas. Dark charcoal background, soft ethereal glow. Minimalist and contemplative.
```

## Content Types

### Announcement Posts
- Lead with the news
- Explain the impact
- Include link or next step

### Insight Posts
- Share one specific learning
- Explain the context briefly
- Make it actionable

### Question Posts
- Ask a genuine question
- Provide your take first
- Keep it focused on one topic

## Quality Checklist

Before finishing:
- [ ] Post saved to `linkedin/<slug>/post.md` or `tweets/<slug>/thread.md`
- [ ] Image generated alongside the post
- [ ] First line hooks attention
- [ ] Content fits platform limits
- [ ] Tone matches platform norms
- [ ] Has clear CTA or question
- [ ] Hashtags are relevant (not generic)

```
See all 185 linesThey instruct the agent to call the `researcher` subagent first, write markdown under `blogs/`, `linkedin/`, or `tweets/`, and call `generate_cover` or `generate_social_image` for images.When you later create the agent and specify the skills folder(s), then the frontmatter of the `SKILLS.md` files from those skill folders get loaded this into the system prompt so the agent can use the skill when a task task matches a skill description.

## [​](#build-the-script)Build the script

Create `content_writer.py` in the project root. The following sections belong in one file, in order.
1[](#)

Add toolsThe researcher subagent uses Tavily search.
Blog and social workflows use Gemini image generation.
When creating the agent later, the `load_subagents` function reads `subagents.yaml` and resolves tool names to these decorated functions.
```
import os
from pathlib import Path
from typing import Literal

import yaml
from langchain.tools import tool

EXAMPLE_DIR = Path(__file__).parent

@tool
def web_search(
 query: str,
 max_results: int = 5,
 topic: Literal["general", "news"] = "general",
) -> dict:
 """Search the web for current information.

 Args:
 query: The search query (be specific and detailed)
 max_results: Number of results to return (default: 5)
 topic: "general" for most queries, "news" for current events

 Returns:
 Search results with titles, URLs, and content excerpts.
 """
 try:
 from tavily import TavilyClient

 api_key = os.environ.get("TAVILY_API_KEY")
 if not api_key:
 return {"error": "TAVILY_API_KEY not set"}

 client = TavilyClient(api_key=api_key)
 return client.search(query, max_results=max_results, topic=topic)
 except Exception as e:
 return {"error": f"Search failed: {e}"}

@tool
def generate_cover(prompt: str, slug: str) -> str:
 """Generate a cover image for a blog post.

 Args:
 prompt: Detailed description of the image to generate.
 slug: Blog post slug. Image saves to blogs/<slug>/hero.png
 """
 try:
 from google import genai

 client = genai.Client()
 response = client.models.generate_content(
 model="gemini-2.5-flash-image",
 contents=[prompt],
 )

 for part in response.parts:
 if part.inline_data is not None:
 image = part.as_image()
 output_path = EXAMPLE_DIR / "blogs" / slug / "hero.png"
 output_path.parent.mkdir(parents=True, exist_ok=True)
 image.save(str(output_path))
 return f"Image saved to {output_path}"

 return "No image generated"
 except Exception as e:
 return f"Error: {e}"

@tool
def generate_social_image(prompt: str, platform: str, slug: str) -> str:
 """Generate an image for a social media post.

 Args:
 prompt: Detailed description of the image to generate.
 platform: Either "linkedin" or "tweets"
 slug: Post slug. Image saves to <platform>/<slug>/image.png
 """
 try:
 from google import genai

 client = genai.Client()
 response = client.models.generate_content(
 model="gemini-2.5-flash-image",
 contents=[prompt],
 )

 for part in response.parts:
 if part.inline_data is not None:
 image = part.as_image()
 output_path = EXAMPLE_DIR / platform / slug / "image.png"
 output_path.parent.mkdir(parents=True, exist_ok=True)
 image.save(str(output_path))
 return f"Image saved to {output_path}"

 return "No image generated"
 except Exception as e:
 return f"Error: {e}"

def load_subagents(config_path: Path) -> list:
 """Load subagent definitions from YAML and wire up tools.

 Unlike `memory` and `skills`, deep agents do not load subagents from files by default.
 This helper externalizes configuration so you can edit YAML without changing Python code.
 """
 available_tools = {
 "web_search": web_search,
 }

 with open(config_path) as f:
 config = yaml.safe_load(f)

 subagents = []
 for name, spec in config.items():
 subagent = {
 "name": name,
 "description": spec["description"],
 "system_prompt": spec["system_prompt"],
 }
 if "model" in spec:
 subagent["model"] = spec["model"]
 if "tools" in spec:
 subagent["tools"] = [available_tools[t] for t in spec["tools"]]
 subagents.append(subagent)

 return subagents

```
2[](#)

Create the agentWhen creating the deep agent with [create_deep_agent](https://reference.langchain.com/python/deepagents/graph/create_deep_agent), pass memory paths, the skills directory, image tools, subagents from YAML, and a [FilesystemBackend](/oss/python/deepagents/backends) rooted at the example directory so paths like `./AGENTS.md` and `./skills/` resolve correctly.GoogleOpenAIAnthropicOpenRouterFireworksBasetenOllama
```
from deepagents import create_deep_agent
from deepagents.backends import FilesystemBackend

def create_content_writer():
 """Create a content writer agent configured by filesystem files."""
 return create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 memory=["./AGENTS.md"],
 skills=["./skills/"],
 tools=[generate_cover, generate_social_image],
 subagents=load_subagents(EXAMPLE_DIR / "subagents.yaml"),
 backend=FilesystemBackend(root_dir=EXAMPLE_DIR),
 )

```
3[](#)

Add an entry pointInvoke the agent with a user message to verify the agent is working:
```
import sys

from langchain.messages import HumanMessage

if __name__ == "__main__":
 task = (
 " ".join(sys.argv[1:])
 if len(sys.argv) > 1
 else "Write a blog post about how AI agents are transforming software development"
 )

 agent = create_content_writer()
 result = agent.invoke(
 {"messages": [HumanMessage(content=task)]},
 config={"configurable": {"thread_id": "content-builder-demo"}},
 )

 for msg in result.get("messages", []):
 if hasattr(msg, "content") and msg.content:
 print(msg.content)

```

## [​](#run-the-agent)Run the agent

The filesystem backend can read, write, and delete files under `root_dir`. Run only in a dedicated directory and review generated content before publishing.
From the project directory you can invoke the agent without passing an argument or by passing the prompt as an argument:
DefaultWith prompt
```
python content_writer.py

```

With `LANGSMITH_API_KEY` set, you can inspect runs in [LangSmith](/langsmith/home).

## [​](#output)Output

On success, generated artifacts are written under a system temporary directory (on macOS and Linux, typically under `/tmp/`), not next to your project files.

```
blogs/
└── prompt-engineering/
 ├── post.md
 └── hero.png
research/
└── prompt-engineering.md

```

Paths follow the skill instructions in `SKILL.md`.

## [​](#full-code)Full code

Browse the complete [content-builder-agent example](https://github.com/langchain-ai/deepagents/tree/main/examples/content-builder-agent) on GitHub, including the Rich-based streaming UI.

## [​](#next-steps)Next steps

{' ' * (self.list_depth - 1)}- Edit `AGENTS.md` to change brand voice and research requirements

{' ' * (self.list_depth - 1)}- Add skills under `skills/<name>/SKILL.md` for new content types

{' ' * (self.list_depth - 1)}- Add subagents in `subagents.yaml` and register tools in `load_subagents`

{' ' * (self.list_depth - 1)}- Read [Subagents](/oss/python/deepagents/subagents), [Skills](/oss/python/deepagents/skills), and [Customization](/oss/python/deepagents/customization) for deeper configuration

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/content-builder.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a deep research agentPrevious](/oss/python/deepagents/deep-research)[Build a data analysis agent from scratchNext](/oss/python/langchain/deep-agent-from-scratch)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
