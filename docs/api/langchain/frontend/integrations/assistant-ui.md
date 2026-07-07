# assistant-ui

[Frontend](/oss/python/langchain/frontend/overview)[Integrations](/oss/python/langchain/frontend/integrations/overview)

# assistant-ui
Copy page

Headless React AI chat framework with a full runtime layer, bridged to useStreamCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.[assistant-ui](https://www.assistant-ui.com/) is a headless React UI framework for AI chat. It provides a full runtime layer—thread management, message branching, attachment handling—that connects to `useStream` via the `useExternalStoreRuntime` adapter.

Clone and run the [full assistant-ui example](https://github.com/langchain-ai/langgraphjs/tree/main/examples/assistant-ui-claude) to see a Claude-style chat interface wired to a LangChain agent with `useExternalStoreRuntime`.

## [​](#how-it-works)How it works

{' ' * (self.list_depth - 1)}- **Stream with `useStream`** — connect to your agent and get reactive messages, loading state, and submit/cancel callbacks

{' ' * (self.list_depth - 1)}- **Adapt with `useExternalStoreRuntime`** — bridge `stream.messages` into assistant-ui’s runtime format by converting `BaseMessage[]` to `ThreadMessageLike[]`

{' ' * (self.list_depth - 1)}- **Provide the runtime** — wrap your UI in `AssistantRuntimeProvider` and render any assistant-ui thread component

## [​](#installation)Installation

```
bun add @assistant-ui/react @assistant-ui/react-markdown

```

## [​](#wiring-usestream)Wiring useStream

The `useExternalStoreRuntime` adapter bridges `stream.messages` into the assistant-ui runtime. Pass it to `AssistantRuntimeProvider` and render any thread component:

```
import { useCallback, useMemo } from "react";
import {
 AssistantRuntimeProvider,
 useExternalStoreRuntime,
 type AppendMessage,
 type ThreadMessageLike,
} from "@assistant-ui/react";
import { useStream } from "@langchain/react";
import { Thread } from "@assistant-ui/react";

export function Chat() {
 const stream = useStream({
 apiUrl: "http://localhost:2024",
 assistantId: "agent",
 });

 const onNew = useCallback(
 async (message: AppendMessage) => {
 const text = message.content
 .filter((c) => c.type === "text")
 .map((c) => c.text)
 .join("");
 await stream.submit({ messages: [{ type: "human", content: text }] });
 },
 [stream],
 );

 // Convert LangChain messages to assistant-ui's ThreadMessageLike format
 const messages = useMemo(
 () => toThreadMessages(stream.messages),
 [stream.messages],
 );

 const runtime = useExternalStoreRuntime<ThreadMessageLike>({
 messages,
 onNew,
 onCancel: () => stream.stop(),
 convertMessage: (m) => m,
 });

 return (
 <AssistantRuntimeProvider runtime={runtime}>
 <Thread />
 </AssistantRuntimeProvider>
 );
}

```

### [​](#converting-messages)Converting messages

`toThreadMessages` maps LangChain `BaseMessage[]` to the `ThreadMessageLike[]` format assistant-ui expects. Handle each message type — human, AI, and tool — and convert content blocks, tool calls, and reasoning tokens:

```
import { AIMessage, HumanMessage, ToolMessage } from "@langchain/core/messages";
import type { ThreadMessageLike } from "@assistant-ui/react";

export function toThreadMessages(messages: BaseMessage[]): ThreadMessageLike[] {
 const result: ThreadMessageLike[] = [];

 for (const msg of messages) {
 if (HumanMessage.isInstance(msg)) {
 result.push({
 role: "user",
 content: [{ type: "text", text: getTextContent(msg.content) }],
 });
 } else if (AIMessage.isInstance(msg)) {
 const parts: ThreadMessageLike["content"] = [];

 // Reasoning tokens
 const reasoning = getReasoningText(msg);
 if (reasoning) parts.push({ type: "reasoning", reasoning });

 // Tool calls
 for (const tc of msg.tool_calls ?? []) {
 parts.push({
 type: "tool-call",
 toolCallId: tc.id ?? "",
 toolName: tc.name,
 args: tc.args,
 });
 }

 // Text response
 const text = getTextContent(msg.content);
 if (text) parts.push({ type: "text", text });

 result.push({ role: "assistant", content: parts });
 } else if (ToolMessage.isInstance(msg)) {
 // Attach tool results to the preceding assistant message
 const last = result[result.length - 1];
 if (last?.role === "assistant") {
 for (const part of last.content) {
 if (
 part.type === "tool-call" &&
 part.toolCallId === msg.tool_call_id
 ) {
 (part as { result?: string }).result = getTextContent(msg.content);
 }
 }
 }
 }
 }

 return result;
}

```
See all 52 lines

## [​](#customising-the-thread-ui)Customising the thread UI

`<Thread />` ships a complete default thread UI including message list, composer, and scroll management. Customise individual parts by overriding component slots:

```
import { Thread, ThreadMessages, Composer } from "@assistant-ui/react";

function CustomThread() {
 return (
 <Thread.Root>
 <ThreadMessages
 components={{
 UserMessage: MyUserMessage,
 AssistantMessage: MyAssistantMessage,
 ToolFallback: MyToolCard,
 }}
 />
 <Composer />
 </Thread.Root>
 );
}

```

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Memoise message conversion:** wrap `toThreadMessages(stream.messages)` in `useMemo` to avoid re-running the conversion on every render

{' ' * (self.list_depth - 1)}- **Handle attachments:** use `CompositeAttachmentAdapter` with `SimpleImageAttachmentAdapter` for image uploads; extend with custom adapters for files

{' ' * (self.list_depth - 1)}- **Use branching:** assistant-ui has built-in message branching support via `MessageBranch`; edit a message to regenerate from that point

{' ' * (self.list_depth - 1)}- **Thread persistence:** `useStream` with `fetchStateHistory: true` and `reconnectOnMount: true` gives assistant-ui access to the full thread history on page load

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/integrations/assistant-ui.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[AI ElementsPrevious](/oss/python/langchain/frontend/integrations/ai-elements)[OpenUINext](/oss/python/langchain/frontend/integrations/openui)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
