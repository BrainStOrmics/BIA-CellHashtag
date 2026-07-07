# Reasoning tokens

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Reasoning tokens
Copy page

Display model thinking and reasoning processes in collapsible blocksCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Reasoning tokens expose the internal thought process of advanced models like OpenAI’s o1/o3 and Anthropic’s Claude with extended thinking. These models produce structured content blocks that separate reasoning from the final answer, letting you build UIs that show *how* the model arrived at its response.

## [​](#what-are-reasoning-tokens)What are reasoning tokens?

When models with reasoning capabilities process a prompt, they generate two distinct types of content:

{' ' * (self.list_depth - 1)}- **Reasoning blocks**: the model’s internal chain-of-thought, problem decomposition, and step-by-step analysis

{' ' * (self.list_depth - 1)}- **Text blocks**: the final, polished response presented to the user

These are delivered as typed content blocks within an `AIMessage`, accessible via the `contentBlocks` property:

```
// Reasoning block
{ type: "reasoning", reasoning: "Let me think about this step by step..." }

// Text block
{ type: "text", text: "The answer is 42." }

```

Not all models produce reasoning tokens. This pattern applies specifically to models that support extended thinking or chain-of-thought output. Standard chat models return only text blocks.

## [​](#use-cases)Use cases

{' ' * (self.list_depth - 1)}- **Transparency**: show users the model’s reasoning process to build trust in its answers

{' ' * (self.list_depth - 1)}- **Debugging**: inspect the model’s thought process to identify where it goes wrong

{' ' * (self.list_depth - 1)}- **Educational tools**: teach students problem-solving by revealing how an AI approaches questions

{' ' * (self.list_depth - 1)}- **Decision support**: let domain experts validate the reasoning behind recommendations

{' ' * (self.list_depth - 1)}- **Quality assurance**: audit reasoning chains for compliance in regulated industries

## [​](#extracting-reasoning-and-text-blocks)Extracting reasoning and text blocks

The `contentBlocks` array on an `AIMessage` contains all blocks in the order they were generated. Filter them by `type` to separate reasoning from text:

```
import { AIMessage } from "@langchain/core/messages";

function extractBlocks(msg: AIMessage) {
 const reasoningBlocks = msg.contentBlocks
 .filter((b) => b.type === "reasoning")
 .map((b) => b.reasoning);

 const textBlocks = msg.contentBlocks
 .filter((b) => b.type === "text")
 .map((b) => b.text);

 return {
 reasoning: reasoningBlocks.join(""),
 text: textBlocks.join(""),
 };
}

```

A single message may contain multiple reasoning blocks (e.g., if the model pauses its reasoning, produces partial text, then reasons further). Joining them gives you the complete thought process.

## [​](#accessing-messages-from-usestream)Accessing messages from `useStream`

Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface AgentState {
 messages: BaseMessage[];
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";
import { AIMessage, HumanMessage } from "@langchain/core/messages";

function Chat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "reasoning",
 });

 return (
 <div className="messages">
 {stream.messages.map((msg, i) => {
 if (HumanMessage.isInstance(msg)) {
 return <HumanBubble key={i} text={msg.content} />;
 }
 if (AIMessage.isInstance(msg)) {
 return (
 <AIResponse
 key={i}
 message={msg}
 isStreaming={stream.isLoading && i === stream.messages.length - 1}
 />
 );
 }
 return null;
 })}
 </div>
 );
}

```

## [​](#building-a-thinkingbubble-component)Building a ThinkingBubble component

The `ThinkingBubble` presents reasoning tokens in a visually distinct, collapsible container. Users can expand it to see the full thought process or collapse it to focus on the final answer.

```
import { useState } from "react";

function ThinkingBubble({
 reasoning,
 isStreaming,
}: {
 reasoning: string;
 isStreaming: boolean;
}) {
 const [isExpanded, setIsExpanded] = useState(false);

 const charCount = reasoning.length;
 const previewLength = 120;
 const preview =
 reasoning.length > previewLength
 ? reasoning.slice(0, previewLength) + "..."
 : reasoning;

 return (
 <div className="thinking-bubble">
 <button
 className="thinking-header"
 onClick={() => setIsExpanded(!isExpanded)}
 >
 <span className="thinking-icon">
 {isStreaming ? (
 <span className="thinking-spinner" />
 ) : (
 "💭"
 )}
 </span>
 <span className="thinking-label">
 {isStreaming ? "Thinking..." : `Thought process (${charCount} chars)`}
 </span>
 <span className={`chevron ${isExpanded ? "expanded" : ""}`}>▶</span>
 </button>

 {isExpanded && (
 <div className="thinking-content">
 <pre>{reasoning}</pre>
 </div>
 )}

 {!isExpanded && !isStreaming && (
 <div className="thinking-preview">{preview}</div>
 )}
 </div>
 );
}

```

### [​](#styling-the-thinkingbubble)Styling the ThinkingBubble

Differentiate reasoning blocks from regular messages with a distinct visual treatment:

```
.thinking-bubble {
 background-color: #f8f5ff;
 border: 1px solid #e2d9f3;
 border-radius: 8px;
 padding: 12px;
 margin: 8px 0;
 font-size: 0.9em;
}

.thinking-header {
 display: flex;
 align-items: center;
 gap: 8px;
 cursor: pointer;
 background: none;
 border: none;
 width: 100%;
 text-align: left;
 color: #6b21a8;
 font-weight: 500;
}

.thinking-content {
 margin-top: 8px;
 padding-top: 8px;
 border-top: 1px solid #e2d9f3;
 white-space: pre-wrap;
 color: #4a4a4a;
 line-height: 1.5;
}

.thinking-preview {
 margin-top: 4px;
 color: #9ca3af;
 font-style: italic;
 font-size: 0.85em;
}

.chevron {
 margin-left: auto;
 transition: transform 0.2s;
}

.chevron.expanded {
 transform: rotate(90deg);
}

```

## [​](#streaming-indicator-for-reasoning)Streaming indicator for reasoning

While the model is still generating reasoning tokens, show an animated indicator to communicate that thinking is in progress:

```
.thinking-spinner {
 display: inline-block;
 width: 16px;
 height: 16px;
 border: 2px solid #e2d9f3;
 border-top-color: #6b21a8;
 border-radius: 50%;
 animation: spin 0.8s linear infinite;
}

@keyframes spin {
 to { transform: rotate(360deg); }
}

```

During streaming, keep the ThinkingBubble collapsed by default and show only the spinner. Expanding mid-stream can cause layout jitter as new tokens arrive. Let users expand after the reasoning phase completes.

## [​](#rendering-the-complete-ai-response)Rendering the complete AI response

Combine the `ThinkingBubble` and a standard text bubble into a single `AIResponse` component:

```
function AIResponse({
 message,
 isStreaming,
}: {
 message: AIMessage;
 isStreaming: boolean;
}) {
 const reasoningBlocks = message.contentBlocks
 .filter((b) => b.type === "reasoning")
 .map((b) => b.reasoning)
 .join("");

 const textBlocks = message.contentBlocks
 .filter((b) => b.type === "text")
 .map((b) => b.text)
 .join("");

 const hasReasoning = reasoningBlocks.length > 0;
 const hasText = textBlocks.length > 0;

 const isReasoningPhase = isStreaming && !hasText;
 const isTextPhase = isStreaming && hasText;

 return (
 <div className="ai-response">
 {hasReasoning && (
 <ThinkingBubble
 reasoning={reasoningBlocks}
 isStreaming={isReasoningPhase}
 />
 )}
 {hasText && (
 <div className="ai-text-bubble">
 <p>{textBlocks}</p>
 {isTextPhase && <span className="cursor-blink">▊</span>}
 </div>
 )}
 </div>
 );
}

```

## [​](#handling-edge-cases)Handling edge cases

### [​](#messages-without-reasoning)Messages without reasoning

Not every AI message will contain reasoning blocks. When `contentBlocks` has only text blocks, render a standard message bubble without the ThinkingBubble.

### [​](#empty-reasoning-blocks)Empty reasoning blocks

Some models produce empty reasoning blocks as placeholders. Filter these out:

```
const meaningfulReasoning = message.contentBlocks
 .filter((b) => b.type === "reasoning" && b.reasoning.trim().length > 0);

```

### [​](#multiple-reasoning-text-cycles)Multiple reasoning-text cycles

A single message can alternate between reasoning and text blocks. If you need to preserve this interleaving, iterate `contentBlocks` in order rather than grouping by type:

```
message.contentBlocks.forEach((block) => {
 if (block.type === "reasoning") {
 // Render ThinkingBubble
 } else if (block.type === "text") {
 // Render text paragraph
 }
});

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Default to collapsed**: show reasoning on demand, not by default

{' ' * (self.list_depth - 1)}- **Show character count**: gives users a quick sense of how much thinking went into the response

{' ' * (self.list_depth - 1)}- **Differentiate visually**: use distinct colors, borders, or backgrounds so reasoning is never confused with the actual answer

{' ' * (self.list_depth - 1)}- **Animate transitions**: smooth expand/collapse animations improve perceived quality

{' ' * (self.list_depth - 1)}- **Consider accessibility**: use proper ARIA attributes (`aria-expanded`, `aria-controls`) on the toggle button

{' ' * (self.list_depth - 1)}- **Truncate in previews**: show a short preview of the reasoning when collapsed so users can decide whether to expand

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/reasoning-tokens.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Branching chatPrevious](/oss/python/langchain/frontend/branching-chat)[Structured outputNext](/oss/python/langchain/frontend/structured-output)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
