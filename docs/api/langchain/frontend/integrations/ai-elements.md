# AI Elements

[Frontend](/oss/python/langchain/frontend/overview)[Integrations](/oss/python/langchain/frontend/integrations/overview)

# AI Elements
Copy page

Composable shadcn/ui-based components for AI chat interfaces with useStreamCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.[AI Elements](https://elements.ai-sdk.dev/) is a composable, shadcn/ui-based component library purpose-built for AI chat interfaces. Components like `Conversation`, `Message`, `Tool`, `Reasoning`, and `PromptInput` are designed to drop directly into any React project and wire to `stream.messages` with minimal glue code.

Clone and run the [full AI Elements example](https://github.com/langchain-ai/langgraphjs/tree/main/examples/ai-elements) to see tool call rendering, reasoning display, streaming messages, and more in a working project.

## [​](#how-it-works)How it works

{' ' * (self.list_depth - 1)}- **Install components as source files:** AI Elements ships via a CLI that adds components directly to your project (shadcn/ui registry style)

{' ' * (self.list_depth - 1)}- **Map messages to components:** iterate `stream.messages`, render `HumanMessage` instances as user bubbles and `AIMessage` instances as assistant responses

{' ' * (self.list_depth - 1)}- **Compose richer UIs:** wrap tool calls in `<Tool>`, reasoning in `<Reasoning>`, and everything in `<Conversation>` for scroll management

## [​](#installation)Installation

Install AI Elements components via the CLI. They’re added as editable source files into your project:

```
npm install @langchain/react @ai-elements/react
npx ai-elements@latest add conversation message prompt-input tool reasoning suggestion

```

## [​](#wiring-usestream)Wiring useStream

Render AI Elements components directly from `stream.messages`. Each LangChain `BaseMessage` maps to a component:

```
import { useStream } from "@langchain/react";
import { HumanMessage, AIMessage } from "@langchain/core/messages";

import {
 Conversation,
 ConversationContent,
 ConversationScrollButton,
} from "@/components/ai-elements/conversation";
import {
 Message,
 MessageContent,
 MessageResponse,
} from "@/components/ai-elements/message";
import {
 Tool,
 ToolHeader,
 ToolContent,
 ToolInput,
 ToolOutput,
} from "@/components/ai-elements/tool";
import {
 Reasoning,
 ReasoningTrigger,
 ReasoningContent,
} from "@/components/ai-elements/reasoning";
import {
 PromptInput,
 PromptInputBody,
 PromptInputTextarea,
 PromptInputFooter,
 PromptInputSubmit,
} from "@/components/ai-elements/prompt-input";

export function Chat() {
 const stream = useStream({
 apiUrl: "http://localhost:2024",
 assistantId: "agent",
 });

 return (
 <div className="flex flex-col h-dvh">
 <Conversation className="flex-1">
 <ConversationContent>
 {stream.messages.map((msg, i) => {
 if (HumanMessage.isInstance(msg)) {
 return (
 <Message key={i} from="user">
 <MessageContent>{msg.content as string}</MessageContent>
 </Message>
 );
 }
 if (AIMessage.isInstance(msg)) {
 return (
 <div key={i}>
 {/* Reasoning block (shows when model emits thinking tokens) */}
 <Reasoning>
 <ReasoningTrigger />
 <ReasoningContent>{getReasoningText(msg)}</ReasoningContent>
 </Reasoning>

 {/* Inline tool calls with input/output display */}
 {getToolCalls(msg).map((tc) => (
 <Tool key={tc.id} defaultOpen>
 <ToolHeader type={`tool-${tc.name}`} state={tc.state} />
 <ToolContent>
 <ToolInput input={tc.args} />
 {tc.output && (
 <ToolOutput output={tc.output} errorText={undefined} />
 )}
 </ToolContent>
 </Tool>
 ))}

 {/* Streamed text response */}
 <Message from="assistant">
 <MessageContent>
 <MessageResponse>{getTextContent(msg)}</MessageResponse>
 </MessageContent>
 </Message>
 </div>
 );
 }
 })}
 </ConversationContent>
 <ConversationScrollButton />
 </Conversation>

 <PromptInput
 onSubmit={({ text }) =>
 stream.submit({ messages: [{ type: "human", content: text }] })
 }
 >
 <PromptInputBody>
 <PromptInputTextarea placeholder="Ask me something..." />
 </PromptInputBody>
 <PromptInputFooter>
 <PromptInputSubmit
 status={stream.isLoading ? "streaming" : "ready"}
 />
 </PromptInputFooter>
 </PromptInput>
 </div>
 );
}

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Edit source files freely:** components ship in your project, not as an external package dependency, so you can change anything without forking

{' ' * (self.list_depth - 1)}- **Use `MessageResponse` for streaming:** it handles streamed partial tokens correctly; avoid rendering raw `msg.content` directly during streaming

{' ' * (self.list_depth - 1)}- **Wrap in `Conversation`:** the `Conversation` component manages scroll behaviour so new messages auto-scroll into view

{' ' * (self.list_depth - 1)}- **Gate on `isInstance`:** use `HumanMessage.isInstance(msg)` and `AIMessage.isInstance(msg)` rather than checking `msg.getType()` for proper TypeScript narrowing

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/integrations/ai-elements.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[CopilotKitPrevious](/oss/python/langchain/frontend/integrations/copilotkit)[assistant-uiNext](/oss/python/langchain/frontend/integrations/assistant-ui)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
